import os
import re
import subprocess
import logging
from typing import Optional, Dict, List, Any


class JobSchedulerError(Exception):
    """Custom exception for scheduler-related errors."""
    pass


class JobScheduler:
    """
    A compact unified interface for job submission, status checking, and cancellation
    across multiple job schedulers (SLURM, PBS, LSF, SGE, HTCondor).
    """

    SCHEDULER_CONFIG = {
        "slurm": {
            "submit": ["sbatch"],
            "status": ["squeue", "--job"],
            "cancel": ["scancel"],
            "id_regex": r"Submitted batch job (\d+)",
            "id_format": r"^\d+$",
            "test_cmd": "sbatch"
        },
        "pbs": {
            "submit": ["qsub"],
            "status": ["qstat"],
            "cancel": ["qdel"],
            "id_regex": r"(\d+(?:\.\S+)?)",
            "id_format": r"^\d+(?:\.\S+)?$",
            "test_cmd": "qsub"
        },
        "torque": {
            "submit": ["qsub"],
            "status": ["qstat"],
            "cancel": ["qdel"],
            "id_regex": r"(\d+(?:\.\S+)?)",
            "id_format": r"^\d+(?:\.\S+)?$",
            "test_cmd": "qsub"
        },
        "lsf": {
            "submit": ["bsub"],
            "status": ["bjobs"],
            "cancel": ["bkill"],
            "id_regex": r"Job <(\d+)> is submitted",
            "id_format": r"^\d+$",
            "test_cmd": "bsub",
            "submit_method": "stdin"  # Special flag for LSF
        },
        "sge": {
            "submit": ["qsub"],
            "status": ["qstat", "-j"],
            "cancel": ["qdel"],
            "id_regex": r"Your job (\d+)",
            "id_format": r"^\d+$",
            "test_cmd": "qsub"
        },
        "htcondor": {
            "submit": ["condor_submit"],
            "status": ["condor_q"],
            "cancel": ["condor_rm"],
            "id_regex": r"submitted to cluster (\d+)",
            "id_format": r"^\d+$",
            "test_cmd": "condor_submit"
        },
        "condor": {  # Alias for htcondor
            "submit": ["condor_submit"],
            "status": ["condor_q"],
            "cancel": ["condor_rm"],
            "id_regex": r"submitted to cluster (\d+)",
            "id_format": r"^\d+$",
            "test_cmd": "condor_submit"
        }
    }

    def __init__(self, scheduler: str, logger: Optional[logging.Logger] = None, timeout: int = 30):
        """
        Initialize the JobScheduler.

        Args:
            scheduler: Name of the scheduler (slurm, pbs, lsf, sge, htcondor)
            logger: Optional logger instance
            timeout: Timeout for subprocess calls in seconds
        """
        self.scheduler_name = scheduler.lower()
        self.logger = logger or logging.getLogger(__name__)
        self.timeout = timeout

        if self.scheduler_name not in self.SCHEDULER_CONFIG:
            raise JobSchedulerError(
                f"Scheduler '{scheduler}' is not supported. "
                f"Supported schedulers: {list(self.SCHEDULER_CONFIG.keys())}"
            )

        self.config = self.SCHEDULER_CONFIG[self.scheduler_name]
        self._validate_scheduler_availability()

    def _validate_scheduler_availability(self) -> None:
        """Check if the scheduler commands are available in the system."""
        test_cmd = self.config["test_cmd"]
        try:
            subprocess.run([test_cmd, "--version"],
                         capture_output=True,
                         timeout=5,
                         check=False)
        except FileNotFoundError:
            raise JobSchedulerError(
                f"Scheduler command '{test_cmd}' not found. "
                f"Please ensure {self.scheduler_name} is installed and in PATH."
            )
        except subprocess.TimeoutExpired:
            # Some schedulers might not support --version, but command exists
            pass

    def _run_command(self, cmd: List[str], input_data: Optional[str] = None) -> subprocess.CompletedProcess:
        """
        Run a command safely with proper error handling.

        Args:
            cmd: Command to execute as a list
            input_data: Optional input data for stdin

        Returns:
            CompletedProcess result

        Raises:
            JobSchedulerError: If command execution fails
        """
        try:
            self.logger.debug(f"Executing command: {' '.join(cmd)}")

            result = subprocess.run(
                cmd,
                input=input_data,
                capture_output=True,
                text=True,
                timeout=self.timeout,
                check=True
            )

            return result

        except subprocess.CalledProcessError as e:
            raise JobSchedulerError(
                f"Command failed with return code {e.returncode}: {e.stderr.strip()}"
            )
        except subprocess.TimeoutExpired:
            raise JobSchedulerError(f"Command timed out after {self.timeout} seconds")
        except FileNotFoundError:
            raise JobSchedulerError(f"Command '{cmd[0]}' not found")
        except Exception as e:
            raise JobSchedulerError(f"Unexpected error: {str(e)}")

    def _extract_job_id(self, output: str) -> Optional[str]:
        """Extract job ID from scheduler output using regex."""
        regex = self.config["id_regex"]
        match = re.search(regex, output)
        return match.group(1) if match else None

    def _validate_job_id(self, job_id: str) -> bool:
        """Validate job ID format for the current scheduler."""
        pattern = self.config["id_format"]
        return bool(re.match(pattern, job_id))

    def _build_command(self, operation: str, argument: Optional[str] = None) -> List[str]:
        """
        Build command for the specified operation.

        Args:
            operation: Operation type (submit, status, cancel)
            argument: Additional argument (script path or job ID)

        Returns:
            Complete command as list
        """
        cmd = self.config[operation].copy()

        if argument:
            cmd.append(argument)

        return cmd

    def submit_job(self, script_path: str, **kwargs) -> str:
        """
        Submit a job script and return the job ID.

        Args:
            script_path: Path to the job script
            **kwargs: Additional scheduler-specific options (future extension)

        Returns:
            Job ID as string

        Raises:
            JobSchedulerError: If submission fails
        """
        if not script_path or not isinstance(script_path, str):
            raise JobSchedulerError("Script path must be a non-empty string")

        script_path = os.path.abspath(script_path)
        if not os.path.exists(script_path):
            raise JobSchedulerError(f"Job script '{script_path}' not found")

        if not os.access(script_path, os.R_OK):
            raise JobSchedulerError(f"Job script '{script_path}' is not readable")

        # Handle special case for LSF (requires stdin input)
        if self.config.get("submit_method") == "stdin":
            cmd = self._build_command("submit")
            with open(script_path, 'r') as f:
                script_content = f.read()
            result = self._run_command(cmd, input_data=script_content)
        else:
            cmd = self._build_command("submit", script_path)
            result = self._run_command(cmd)

        job_id = self._extract_job_id(result.stdout)
        if not job_id:
            raise JobSchedulerError(
                f"Could not extract job ID from output: {result.stdout.strip()}"
            )

        self.logger.info(f"Job submitted successfully with ID: {job_id}")
        return job_id

    def check_status(self, job_id: str) -> str:
        """
        Check the status of a job.

        Args:
            job_id: Job ID to check

        Returns:
            Status information as string

        Raises:
            JobSchedulerError: If status check fails
        """
        if not job_id:
            raise JobSchedulerError("Job ID cannot be empty")

        if not self._validate_job_id(job_id):
            raise JobSchedulerError(f"Invalid job ID format: {job_id}")

        cmd = self._build_command("status", job_id)

        try:
            result = self._run_command(cmd)
            return result.stdout.strip()
        except JobSchedulerError as e:
            # Some schedulers return non-zero for completed jobs
            if "return code" in str(e):
                return f"Job {job_id} may be completed or not found"
            raise

    def cancel_job(self, job_id: str) -> bool:
        """
        Cancel a running job.

        Args:
            job_id: Job ID to cancel

        Returns:
            True if cancellation successful, False otherwise

        Raises:
            JobSchedulerError: If cancellation fails
        """
        if not job_id:
            raise JobSchedulerError("Job ID cannot be empty")

        if not self._validate_job_id(job_id):
            raise JobSchedulerError(f"Invalid job ID format: {job_id}")

        cmd = self._build_command("cancel", job_id)

        try:
            self._run_command(cmd)
            self.logger.info(f"Job {job_id} cancelled successfully")
            return True
        except JobSchedulerError as e:
            self.logger.error(f"Failed to cancel job {job_id}: {str(e)}")
            return False

    def get_supported_schedulers(self) -> List[str]:
        """Return list of supported scheduler names."""
        return list(self.SCHEDULER_CONFIG.keys())

    def get_scheduler_info(self) -> Dict[str, Any]:
        """Return configuration info for current scheduler."""
        return {
            "name": self.scheduler_name,
            "commands": {
                "submit": self.config["submit"],
                "status": self.config["status"],
                "cancel": self.config["cancel"]
            },
            "id_format": self.config["id_format"]
        }

    def __repr__(self) -> str:
        return f"JobScheduler(scheduler='{self.scheduler_name}')"


# Example usage and testing
if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    try:
        # Initialize scheduler
        scheduler = JobScheduler("slurm")

        # Print scheduler info
        print(f"Scheduler info: {scheduler.get_scheduler_info()}")
        print(f"Supported schedulers: {scheduler.get_supported_schedulers()}")

        # Example usage (commented out to avoid actual execution)
        # job_id = scheduler.submit_job("my_script.sh")
        # print(f"Submitted job: {job_id}")

        # status = scheduler.check_status(job_id)
        # print(f"Job status: {status}")

        # success = scheduler.cancel_job(job_id)
        # print(f"Cancellation successful: {success}")

    except JobSchedulerError as e:
        print(f"Error: {e}")
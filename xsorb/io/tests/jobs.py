import os
import re
import subprocess
import logging
from typing import Optional, Dict, List, Union
from abc import ABC, abstractmethod


class JobSchedulerError(Exception):
    """Custom exception for scheduler-related errors."""
    pass


class BaseScheduler(ABC):
    """Abstract base class for job schedulers."""

    @abstractmethod
    def get_submit_command(self, script_path: str) -> List[str]:
        pass

    @abstractmethod
    def get_status_command(self, job_id: str) -> List[str]:
        pass

    @abstractmethod
    def get_cancel_command(self, job_id: str) -> List[str]:
        pass

    @abstractmethod
    def extract_job_id(self, output: str) -> Optional[str]:
        pass

    @abstractmethod
    def validate_job_id(self, job_id: str) -> bool:
        pass


class SlurmScheduler(BaseScheduler):
    """SLURM scheduler implementation."""

    def get_submit_command(self, script_path: str) -> List[str]:
        return ["sbatch", script_path]

    def get_status_command(self, job_id: str) -> List[str]:
        return ["squeue", "--job", job_id, "--format=%T,%R"]

    def get_cancel_command(self, job_id: str) -> List[str]:
        return ["scancel", job_id]

    def extract_job_id(self, output: str) -> Optional[str]:
        match = re.search(r"Submitted batch job (\d+)", output)
        return match.group(1) if match else None

    def validate_job_id(self, job_id: str) -> bool:
        return job_id.isdigit()


class PbsScheduler(BaseScheduler):
    """PBS/Torque scheduler implementation."""

    def get_submit_command(self, script_path: str) -> List[str]:
        return ["qsub", script_path]

    def get_status_command(self, job_id: str) -> List[str]:
        return ["qstat", job_id]

    def get_cancel_command(self, job_id: str) -> List[str]:
        return ["qdel", job_id]

    def extract_job_id(self, output: str) -> Optional[str]:
        # PBS typically returns: "123.hostname" or just "123"
        match = re.search(r"(\d+(?:\.\S+)?)", output.strip())
        return match.group(1) if match else None

    def validate_job_id(self, job_id: str) -> bool:
        return bool(re.match(r"\d+(?:\.\S+)?$", job_id))


class LsfScheduler(BaseScheduler):
    """LSF scheduler implementation."""

    def get_submit_command(self, script_path: str) -> List[str]:
        return ["bsub", "<", script_path]

    def get_status_command(self, job_id: str) -> List[str]:
        return ["bjobs", job_id]

    def get_cancel_command(self, job_id: str) -> List[str]:
        return ["bkill", job_id]

    def extract_job_id(self, output: str) -> Optional[str]:
        match = re.search(r"Job <(\d+)> is submitted", output)
        return match.group(1) if match else None

    def validate_job_id(self, job_id: str) -> bool:
        return job_id.isdigit()


class SgeScheduler(BaseScheduler):
    """SGE (Sun Grid Engine) scheduler implementation."""

    def get_submit_command(self, script_path: str) -> List[str]:
        return ["qsub", script_path]

    def get_status_command(self, job_id: str) -> List[str]:
        return ["qstat", "-j", job_id]

    def get_cancel_command(self, job_id: str) -> List[str]:
        return ["qdel", job_id]

    def extract_job_id(self, output: str) -> Optional[str]:
        match = re.search(r"Your job (\d+)", output)
        return match.group(1) if match else None

    def validate_job_id(self, job_id: str) -> bool:
        return job_id.isdigit()


class HtCondorScheduler(BaseScheduler):
    """HTCondor scheduler implementation."""

    def get_submit_command(self, script_path: str) -> List[str]:
        return ["condor_submit", script_path]

    def get_status_command(self, job_id: str) -> List[str]:
        return ["condor_q", job_id]

    def get_cancel_command(self, job_id: str) -> List[str]:
        return ["condor_rm", job_id]

    def extract_job_id(self, output: str) -> Optional[str]:
        # HTCondor: "1 job(s) submitted to cluster 123."
        match = re.search(r"submitted to cluster (\d+)", output)
        return match.group(1) if match else None

    def validate_job_id(self, job_id: str) -> bool:
        return job_id.isdigit()


class JobScheduler:
    """
    A unified interface for job submission, status checking, and cancellation
    across multiple job schedulers (SLURM, PBS, LSF, SGE, HTCondor).
    """

    SUPPORTED_SCHEDULERS = {
        "slurm": SlurmScheduler,
        "pbs": PbsScheduler,
        "torque": PbsScheduler,  # Torque uses PBS commands
        "lsf": LsfScheduler,
        "sge": SgeScheduler,
        "htcondor": HtCondorScheduler,
        "condor": HtCondorScheduler,
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

        if self.scheduler_name not in self.SUPPORTED_SCHEDULERS:
            raise JobSchedulerError(
                f"Scheduler '{scheduler}' is not supported. "
                f"Supported schedulers: {list(self.SUPPORTED_SCHEDULERS.keys())}"
            )

        self.scheduler = self.SUPPORTED_SCHEDULERS[self.scheduler_name]()
        self._validate_scheduler_availability()

    def _validate_scheduler_availability(self) -> None:
        """Check if the scheduler commands are available in the system."""
        test_commands = {
            "slurm": "sbatch",
            "pbs": "qsub",
            "torque": "qsub",
            "lsf": "bsub",
            "sge": "qsub",
            "htcondor": "condor_submit",
            "condor": "condor_submit",
        }

        command = test_commands.get(self.scheduler_name)
        if command:
            try:
                subprocess.run([command, "--version"],
                             capture_output=True,
                             timeout=5,
                             check=False)
            except FileNotFoundError:
                raise JobSchedulerError(
                    f"Scheduler command '{command}' not found. "
                    f"Please ensure {self.scheduler_name} is installed and in PATH."
                )
            except subprocess.TimeoutExpired:
                # Some schedulers might not support --version, but command exists
                pass

    def _run_command(self, cmd: List[str]) -> subprocess.CompletedProcess:
        """
        Run a command safely with proper error handling.

        Args:
            cmd: Command to execute as a list

        Returns:
            CompletedProcess result

        Raises:
            JobSchedulerError: If command execution fails
        """
        try:
            self.logger.debug(f"Executing command: {' '.join(cmd)}")

            # Handle LSF special case with input redirection
            if self.scheduler_name == "lsf" and "<" in cmd:
                idx = cmd.index("<")
                script_path = cmd[idx + 1]
                cmd = cmd[:idx]

                with open(script_path, 'r') as f:
                    script_content = f.read()

                result = subprocess.run(
                    cmd,
                    input=script_content,
                    capture_output=True,
                    text=True,
                    timeout=self.timeout,
                    check=True
                )
            else:
                result = subprocess.run(
                    cmd,
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

        cmd = self.scheduler.get_submit_command(script_path)
        result = self._run_command(cmd)

        job_id = self.scheduler.extract_job_id(result.stdout)
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

        if not self.scheduler.validate_job_id(job_id):
            raise JobSchedulerError(f"Invalid job ID format: {job_id}")

        cmd = self.scheduler.get_status_command(job_id)

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

        if not self.scheduler.validate_job_id(job_id):
            raise JobSchedulerError(f"Invalid job ID format: {job_id}")

        cmd = self.scheduler.get_cancel_command(job_id)

        try:
            self._run_command(cmd)
            self.logger.info(f"Job {job_id} cancelled successfully")
            return True
        except JobSchedulerError as e:
            self.logger.error(f"Failed to cancel job {job_id}: {str(e)}")
            return False

    def get_supported_schedulers(self) -> List[str]:
        """Return list of supported scheduler names."""
        return list(self.SUPPORTED_SCHEDULERS.keys())

    def __repr__(self) -> str:
        return f"JobScheduler(scheduler='{self.scheduler_name}')"


# Example usage and testing
if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    try:
        # Initialize scheduler
        scheduler = JobScheduler("slurm")

        # Submit a job (example)
        # job_id = scheduler.submit_job("my_script.sh")
        # print(f"Submitted job: {job_id}")

        # Check status
        # status = scheduler.check_status(job_id)
        # print(f"Job status: {status}")

        # Cancel job
        # success = scheduler.cancel_job(job_id)
        # print(f"Cancellation successful: {success}")

        print(f"Supported schedulers: {scheduler.get_supported_schedulers()}")

    except JobSchedulerError as e:
        print(f"Error: {e}")
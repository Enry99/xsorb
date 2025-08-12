"""
Module to handle job submission, status checking, and cancellation
across multiple job schedulers in a unified way.
"""

import os
import re
import subprocess
import logging
from typing import Optional, List

# Update SCHEDULER_CONFIG to include commands for active jobs
SCHEDULER_CONFIG = {
    "slurm": {
        "submit": ["sbatch"],
        "cancel": ["scancel"],
        "id_regex": r"Submitted batch job (\d+)",
        "id_format": r"^\d+$",
        "test_cmd": "sbatch",
        "list_active_jobs": ["squeue", "-h", "-u", "$USER", "-t", "RUNNING,PENDING"]  # Running + queued
    },
    "pbs": {
        "submit": ["qsub"],
        "cancel": ["qdel"],
        "id_regex": r"(\d+(?:\.\S+)?)",
        "id_format": r"^\d+(?:\.\S+)?$",
        "test_cmd": "qsub",
        "list_active_jobs": ["qstat", "-u", "$USER"]  # All active jobs (R + Q states)
    },
    "torque": {
        "submit": ["qsub"],
        "cancel": ["qdel"],
        "id_regex": r"(\d+(?:\.\S+)?)",
        "id_format": r"^\d+(?:\.\S+)?$",
        "test_cmd": "qsub",
        "list_active_jobs": ["qstat", "-u", "$USER"]  # All active jobs (R + Q states)
    },
    "lsf": {
        "submit": ["bsub"],
        "cancel": ["bkill"],
        "id_regex": r"Job <(\d+)> is submitted",
        "id_format": r"^\d+$",
        "test_cmd": "bsub",
        "submit_method": "stdin",
        "list_active_jobs": ["bjobs", "-u", "$USER"]  # Running + pending jobs
    },
    "sge": {
        "submit": ["qsub"],
        "cancel": ["qdel"],
        "id_regex": r"Your job (\d+)",
        "id_format": r"^\d+$",
        "test_cmd": "qsub",
        "list_active_jobs": ["qstat", "-u", "$USER", "-s", "r,qw"]  # Running + queued
    },
    "htcondor": {
        "submit": ["condor_submit"],
        "cancel": ["condor_rm"],
        "id_regex": r"submitted to cluster (\d+)",
        "id_format": r"^\d+$",
        "test_cmd": "condor_submit",
        "list_active_jobs": ["condor_q", "$USER"]  # All active jobs (idle + running)
    },
    "condor": {  # Alias for htcondor
        "submit": ["condor_submit"],
        "cancel": ["condor_rm"],
        "id_regex": r"submitted to cluster (\d+)",
        "id_format": r"^\d+$",
        "test_cmd": "condor_submit",
        "list_active_jobs": ["condor_q", "$USER"]  # All active jobs (idle + running)
    }
}


class JobScheduler:
    """
    A compact unified interface for job submission, status checking, and cancellation
    across multiple job schedulers (SLURM, PBS, LSF, SGE, HTCondor).
    """

    def __init__(self, scheduler: str | dict,
                 logger: Optional[logging.Logger] = None,
                 timeout: int = 30):
        """
        Initialize the JobScheduler.

        Args:
            scheduler: Name of the scheduler (slurm, pbs, lsf, sge, htcondor),
                        or a custom configuration dictionary
            logger: Optional logger instance
            timeout: Timeout for subprocess calls in seconds
        """

        if isinstance(scheduler, dict):
            SCHEDULER_CONFIG["custom"] = scheduler
            scheduler = "custom"

        self.scheduler_name = scheduler.lower()
        self.logger = logger or logging.getLogger(__name__)
        self.timeout = timeout

        if self.scheduler_name not in SCHEDULER_CONFIG:
            raise JobSchedulerError(
                f"Scheduler '{scheduler}' is not supported. "
                f"Supported schedulers: {list(SCHEDULER_CONFIG.keys())}. "
                "Please provide a valid scheduler name or a custom configuration."
            )

        self.config = SCHEDULER_CONFIG[self.scheduler_name]

    def _run_command(self,
                     cmd: List[str],
                     input_data: Optional[str] = None) -> subprocess.CompletedProcess:
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
            ) from e
        except subprocess.TimeoutExpired as e:
            raise JobSchedulerError(f"Command timed out after {self.timeout} seconds") from e
        except FileNotFoundError as e:
            raise JobSchedulerError(f"Command '{cmd[0]}' not found") from e
        except Exception as e:
            raise JobSchedulerError(f"Unexpected error: {str(e)}") from e

    def _extract_job_id(self, output: str) -> Optional[int]:
        """Extract job ID from scheduler output using regex."""
        regex = self.config["id_regex"]
        match = re.search(regex, output)
        return int(match.group(1)) if match else None

    def _validate_job_id(self, job_id: str) -> bool:
        """Validate job ID format for the current scheduler."""
        pattern = self.config["id_format"]
        return bool(re.match(pattern, job_id))

    def _build_command(self,
                       operation: str,
                       argument: Optional[str] = None,
                       extra_args: Optional[List[str]] = None) -> List[str]:
        """
        Build command for the specified operation.

        Args:
            operation: Operation type (submit, cancel, list_jobs, list_all_jobs)
            argument: Additional argument (script path or job ID)
            extra_args: Extra arguments to append to the command

        Returns:
            Complete command as list
        """
        cmd = self.config[operation].copy()

        if argument:
            cmd.append(argument)

        if extra_args:
            cmd.extend(extra_args)

        return cmd

    def submit_job(self, script_path: str, script_args: Optional[List[str]] = None) -> int:
        """
        Submit a job script and return the job ID.

        Args:
            script_path: Path to the job script
            script_args: List of arguments to pass to the job script (e.g., ["input.pwi", "output.pwo"])

        Returns:
            Job ID as an integer

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

        # Validate script_args if provided
        if script_args is not None:
            if not isinstance(script_args, list):
                raise JobSchedulerError("script_args must be a list of strings")
            if not all(isinstance(arg, str) for arg in script_args):
                raise JobSchedulerError("All script_args must be strings")

        # Handle special case for LSF (requires stdin input)
        if self.config.get("submit_method") == "stdin":
            cmd = self._build_command("submit")

            # For LSF, we need to modify the script content to include arguments
            with open(script_path, 'r') as f:
                script_content = f.read()

            # If script_args are provided for LSF, we need to handle them differently
            if script_args:
                # Add arguments as environment variables or modify the script
                # This is LSF-specific behavior
                env_vars = "\n".join([f"export ARG{i+1}={arg}" for i, arg in enumerate(script_args)])
                script_content = f"{env_vars}\n{script_content}"
                self.logger.debug(f"LSF: Added environment variables for script arguments: {script_args}")

            result = self._run_command(cmd, input_data=script_content)
        else:
            # For other schedulers, arguments are passed directly to the command
            cmd = self._build_command("submit", script_path, script_args)
            result = self._run_command(cmd)

        job_id = self._extract_job_id(result.stdout)
        if not job_id:
            raise JobSchedulerError(
                f"Could not extract job ID from output: {result.stdout.strip()}"
            )

        log_msg = f"Submitted job {job_id}"
        if script_args:
            log_msg += f" with arguments: {' '.join(script_args)}"
        self.logger.info(log_msg)
        return job_id

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
            self.logger.info(f"Cancelled job {job_id}")
            return True
        except JobSchedulerError as e:
            self.logger.error(f"Failed to cancel job {job_id}: {str(e)}")
            return False

    def get_active_job_ids(self) -> List[str]:
        """
        Get all active job IDs for the current user (running + queued/pending).
        These are jobs that can be cancelled.

        Returns:
            List of job IDs (strings) for jobs that are running or queued

        Raises:
            JobSchedulerError: If command execution fails
        """
        if "list_active_jobs" not in self.config:
            raise JobSchedulerError(f"Active job listing not supported for {self.scheduler_name}")

        cmd = self.config["list_active_jobs"].copy()

        # Replace $USER placeholder with actual username
        if "$USER" in cmd:
            username = os.getenv("USER") or os.getenv("USERNAME")
            if not username:
                raise JobSchedulerError("Cannot determine current username")
            cmd = [arg.replace("$USER", username) if arg == "$USER" else arg for arg in cmd]

        try:
            result = self._run_command(cmd)
            job_ids = self._extract_active_job_ids(result.stdout)

            self.logger.debug(f"Found {len(job_ids)} active jobs for user")
            return job_ids

        except JobSchedulerError as e:
            self.logger.error(f"Failed to get active jobs: {str(e)}")
            raise

    def _extract_active_job_ids(self, output: str) -> List[str]:
        """
        Extract active job IDs from scheduler output (running + queued).

        Args:
            output: Raw output from list active jobs command

        Returns:
            List of job IDs for active jobs
        """
        job_ids = []
        lines = output.strip().split('\n')

        for line in lines:
            if not line.strip():
                continue

            fields = line.split()
            if not fields:
                continue

            # First field is usually the job ID for most schedulers
            job_id = fields[0]

            # Validate job ID format
            if self._validate_job_id(job_id):
                # For schedulers that include completed jobs, filter only active ones
                if self.scheduler_name in ["pbs", "torque"]:
                    # Check status column (typically field 5) - include R (running) and Q (queued)
                    if len(fields) >= 6 and fields[5] in ["R", "Q", "H", "T", "W", "S"]:
                        job_ids.append(job_id)
                elif self.scheduler_name == "lsf":
                    # Check status column (typically field 2) - include RUN and PEND states
                    if len(fields) >= 3 and fields[2] in ["RUN", "PEND", "PSUSP", "USUSP", "SSUSP"]:
                        job_ids.append(job_id)
                else:
                    # For SLURM, SGE, HTCondor - command already filters active jobs
                    job_ids.append(job_id)

        return job_ids

    def get_supported_schedulers(self) -> List[str]:
        """Return list of supported scheduler names."""
        return list(SCHEDULER_CONFIG.keys())

    def __repr__(self) -> str:
        return f"JobScheduler(scheduler='{self.scheduler_name}')"


class JobSchedulerError(Exception):
    """Custom exception for scheduler-related errors."""
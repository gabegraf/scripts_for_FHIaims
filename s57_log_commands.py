import subprocess, datetime, argparse, sys

def runlog(cmd, log="log.txt", append_output=False):
    with open(log, "a", encoding="utf-8") as f:
        f.write(f"[{datetime.datetime.now():%Y-%m-%d %H:%M:%S}] {cmd}\n")
    if append_output:
        p = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        with open(log, "a", encoding="utf-8") as f:
            f.write(p.stdout or "")
            f.write(p.stderr or "")
        # also print captured output to console
        if p.stdout:
            print(p.stdout, end="")
        if p.stderr:
            print(p.stderr, end="", file=sys.stderr)
        return p
    else:
        return subprocess.run(cmd, shell=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run a shell command and log the command (and optionally output)."
    )
    parser.add_argument(
        "-l", "--log", default="log.txt", help="log file path (default: log.txt)"
    )
    parser.add_argument(
        "-a",
        "--append-output",
        action="store_true",
        help="capture stdout/stderr and append to the log",
    )
    parser.add_argument(
        "-c", "--command", help="command string to run (alternative to positional args)"
    )
    parser.add_argument(
        "cmd", nargs=argparse.REMAINDER, help="command and args (use either this or -c)"
    )

    args = parser.parse_args()

    if args.command:
        cmd_str = args.command
    else:
        if not args.cmd:
            parser.error("no command provided")
        # argparse.REMAINDER keeps the raw remaining tokens; join into a single command string
        cmd_str = " ".join(args.cmd).strip()

    result = runlog(cmd_str, log=args.log, append_output=args.append_output)
    # exit with the subprocess return code if available
    sys.exit(
        result.returncode if isinstance(result, subprocess.CompletedProcess) else 0
    )

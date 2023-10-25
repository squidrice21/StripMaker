#!/usr/bin/env python

"""
Utilities for running async jobs in background.
"""

import asyncio
import locale
import sys
import threading

lock = threading.Lock()
processes = []


async def read_stream(stream):
    lines = []
    while not stream.at_eof():
        data = await stream.readline()
        line = data.decode(locale.getpreferredencoding(False))
        if len(line):
            lines.append(line)
    return lines


async def run_subprocess_in_bg(case, cmd, timeout):
    print("Running {}".format(case))
    # https://stackoverflow.com/questions/45769985/asyncio-create-subprocess-exec-console-window-opening-for-each-call
    # This flag avoids popup window after crashing
    if sys.platform == "win32":
        DETACHED_PROCESS = 0x00000008
        process = await asyncio.create_subprocess_exec(
            *cmd,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE,
            creationflags=DETACHED_PROCESS,
        )
    else:
        process = await asyncio.create_subprocess_exec(
            *cmd,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE
        )

    with lock:
        processes.append(process)

    task_code = asyncio.ensure_future(process.wait())
    task_out = asyncio.ensure_future(read_stream(process.stdout))
    task_err = asyncio.ensure_future(read_stream(process.stderr))

    done, pending = await asyncio.wait([task_code, task_out, task_err], timeout=timeout)
    if pending:
        # timeout
        if process.returncode is None:
            # kill the subprocess, then `await future` will return soon
            try:
                print(f'* Killing {case} after {timeout} s...')
                process.kill()
            except ProcessLookupError:
                pass

    code = await task_code
    out = await task_out
    err = await task_err

    if code == 0:
        print(f'* Success: {case}.')
    else:
        print(f'* Failure: {case} {code}...')

    return code, out, err, case


async def run_with_limit(case, cmd, timeout, semaphore):
    # prevent more than certain number things to run at the same time
    async with semaphore:
        return await run_subprocess_in_bg(case, cmd, timeout)


def async_entrance(main_func):
    if sys.platform == "win32":
        loop = asyncio.ProactorEventLoop()
        asyncio.set_event_loop(loop)
    try:
        asyncio.get_event_loop().run_until_complete(main_func())
    except KeyboardInterrupt:
        with lock:
            for proc in processes:
                try:
                    proc.kill()
                except ProcessLookupError:
                    pass

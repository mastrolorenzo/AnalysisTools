import errno
import os


def makedirs(path):
    """Recursively create a directory without race conditions."""
    try:
        os.makedirs(path)
    except OSError as e:
        if e.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


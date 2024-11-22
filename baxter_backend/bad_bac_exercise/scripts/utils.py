import os

def die_if_not_nonempty_file(filename, prepend="", logspecial=None):
    errmsg = None
    if not filename:
        errmsg = f"{prepend} filename not specified."
        if isinstance(filename, str): errmsg += " (An empty string passed as filename.)"
        elif filename is None: errmsg += " (The filename is None.)"
        if isinstance(filename, list): errmsg += " (An empty list passed as filename.)"
    elif not os.path.exists(filename):
        errmsg = f"{prepend} file {filename} not found."
    elif os.path.isdir(filename):
        errmsg = f"{prepend} {filename} is a directory, not a file."
    elif os.path.getsize(filename) == 0:
        errmsg = f"{prepend} {filename} is empty."

    if errmsg:
        raise FileNotFoundError(errmsg)
    return

def is_nonempty_file(fnm):
    return os.path.exists(fnm) and os.path.getsize(fnm) > 0


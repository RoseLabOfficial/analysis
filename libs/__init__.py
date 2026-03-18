def hyperlink(path: str) -> str:
    return f"\u001b]8;;file://{path}\u001b\\{path}\u001b]8;;\u001b\\"
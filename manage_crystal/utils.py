"""Utilities and shortcuts for manage_crystal."""


def is_number(s):
    """Checks if a string is a number or not."""
    try:
        float(s)
        return True
    except ValueError:
        pass
    try:
        import unicodedata
        unicodedata.numeric(s)
        return True
    except (TypeError, ValueError):
        pass


def parse_and_write(inputfile, outputfile):
    from manage_crystal.file_parser import parse_from_filepath
    from manage_crystal.file_writer import write_to_filepath
    crys = parse_from_filepath(inputfile)
    write_to_filepath(crys, outputfile)


def parse_coord(coord_string):
    """Convert string to number, and exclude brackets such as 0.342(3)."""
    coord = coord_string.split("(")[0]
    return float(coord)

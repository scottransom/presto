import subprocess
import json
import os
import os.path
import platform
import sys


# For simple terminal colors
class bcolors:
    HEADER = "\033[95m"
    OKBLUE = "\033[94m"
    OKCYAN = "\033[96m"
    OKGREEN = "\033[92m"
    WARNING = "\033[93m"
    FAIL = "\033[91m"
    ENDC = "\033[0m"
    BOLD = "\033[1m"
    UNDERLINE = "\033[4m"


# Bail with error
def bail():
    print(
        f"\n{bcolors.FAIL}There seem to be some issues. Please fix before trying to build!{bcolors.ENDC}"
    )
    sys.exit(1)


# Run meson introspect and get the key buildoptions as JSON
result = subprocess.run(
    ["meson", "introspect", "build", "--buildoptions"], stdout=subprocess.PIPE
)

# Process the JSON
vars = json.loads(result.stdout)

# Which variables do we want to see?
grab = ["prefix", "libdir", "bindir"]
ourvars = {}

# Grab the important values
for var in vars:
    if var["name"] in grab:
        ourvars[var["name"]] = var["value"]

# Show the values for prefix and libdir
print("\nPRESTO's meson build directory currently has:")
print(f"  prefix = {ourvars['prefix']}")
print(f"  bindir = {ourvars['bindir']}")
print(f"  libdir = {ourvars['libdir']}\n")
libinstall = os.path.join(ourvars["prefix"], ourvars["libdir"])
bininstall = os.path.join(ourvars["prefix"], ourvars["bindir"])

# Now show the user what they should be setting
print(
    f"Checking for {bcolors.BOLD}PRESTO{bcolors.ENDC} as an environment variable:",
    end="",
)
presto = os.environ.get("PRESTO")
if presto is None:
    print(
        f"\n  {bcolors.WARNING}WARNING:{bcolors.ENDC} PRESTO environment variable is not set!"
    )
    print("  You can set it using something like:")
    print(f"  export PRESTO={os.getcwd()}")
    bail()
else:
    print(f"  {bcolors.OKGREEN}yes{bcolors.ENDC}")
    if os.getcwd() != os.path.realpath(presto):
        print(
            f"\n  {bcolors.WARNING}WARNING:{bcolors.ENDC} PRESTO is set, but not to the current directory!"
        )

print(
    f"Checking for {bcolors.BOLD}TEMPO{bcolors.ENDC} as an environment variable:",
    end="",
)
tempo = os.environ.get("TEMPO")
if tempo is None:
    print(f"  {bcolors.FAIL}no{bcolors.ENDC}")
    print(
        f"\n  {bcolors.WARNING}WARNING:{bcolors.ENDC} TEMPO environment variable is not set!"
    )
    print("  You need to set it to the location of the TEMPO code.")
    bail()
else:
    print(f"  {bcolors.OKGREEN}yes{bcolors.ENDC}")

chk = os.path.join(presto, "bin")
path = os.environ.get("PATH")
path = path.split(":") if path is not None else []
path = [x.rstrip("/") for x in path]
if chk != bininstall:
    print("Is $PRESTO/bin in PATH? (it shouldn't be):", end="")
    if chk in path:
        print(f"  {bcolors.FAIL}yes{bcolors.ENDC}")
        print(
            f"\n  {bcolors.WARNING}WARNING:{bcolors.ENDC} $PRESTO/bin should probably not be in your PATH!"
        )
        bail()
    else:
        print(f"  {bcolors.OKGREEN}no{bcolors.ENDC}")

dlpath = "DYLD_LIBRARY_PATH" if platform.system() == "Darwin" else "LD_LIBRARY_PATH"
chk = os.path.join(presto, "lib")
if chk != libinstall:
    print(f"Is $PRESTO/lib in {dlpath}? (it shouldn't be):", end="")
    if chk in path:
        print(f"  {bcolors.FAIL}yes{bcolors.ENDC}")
        print(
            f"\n  {bcolors.WARNING}WARNING:{bcolors.ENDC} $PRESTO/lib should probably not be in your {dlpath}!"
        )
        bail()
    else:
        print(f"  {bcolors.OKGREEN}no{bcolors.ENDC}")

print(f"Is {bcolors.BOLD}{bininstall}{bcolors.ENDC} in PATH? (it should be):", end="")
if bininstall not in path:
    print(f"  {bcolors.FAIL}no{bcolors.ENDC}")
    print(f"\n  {bcolors.WARNING}WARNING:{bcolors.ENDC} {bininstall} is not in PATH!")
    print(
        f"  That may cause issues with linking if {bininstall} is not a standard binary directory."
    )
    print("  You can ensure that it is used by setting PATH using something like:")
    print(f"  export PATH={bininstall}:$PATH")
    bail()
else:
    print(f"  {bcolors.OKGREEN}yes{bcolors.ENDC}")

print(
    f"Is {bcolors.BOLD}{libinstall}{bcolors.ENDC} in LIBRARY_PATH? (it probably should be)",
    end="",
)
libpath = os.environ.get("LIBRARY_PATH")
libpath = libpath.split(":") if libpath is not None else []
libpath = [x.rstrip("/") for x in libpath]
if libinstall not in libpath:
    print(f"  {bcolors.FAIL}no{bcolors.ENDC}")
    print(
        f"\n  {bcolors.WARNING}WARNING:{bcolors.ENDC} {libinstall} is not in LIBRARY_PATH!"
    )
    print(
        f"  That may cause link issues if {libinstall} is not a standard library directory."
    )
    print(
        "  You can ensure that it is used by setting LIBRARY_PATH using something like:"
    )
    print(f"  export LIBRARY_PATH={libinstall}:$LIBRARY_PATH")
    bail()
else:
    print(f"  {bcolors.OKGREEN}yes{bcolors.ENDC}")

print(
    f"Is {bcolors.BOLD}{libinstall}{bcolors.ENDC} in {dlpath}? (it probably should be)",
    end="",
)
ldlibpath = os.environ.get(dlpath)
ldlibpath = ldlibpath.split(":") if ldlibpath is not None else []
ldlibpath = [x.rstrip("/") for x in ldlibpath]
if libinstall not in ldlibpath:
    print(f"  {bcolors.FAIL}no{bcolors.ENDC}")
    print(
        f"\n  {bcolors.WARNING}WARNING:{bcolors.ENDC} {libinstall} is not in {dlpath}!"
    )
    print(
        f"  That may cause runtime issues if {libinstall} is not a standard library directory."
    )
    print(f"  You can ensure that it is used by setting {dlpath} using something like:")
    print(f"  export {dlpath}={libinstall}:${dlpath}")
    bail()
else:
    print(f"  {bcolors.OKGREEN}yes{bcolors.ENDC}")

print(f"\n{bcolors.OKGREEN}Everything looks good! Let's try to build:{bcolors.ENDC}\n"
    "  meson compile -C build\n  meson install -C build")

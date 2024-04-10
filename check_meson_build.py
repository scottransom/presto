import subprocess
import json
import os
import os.path

# For simple terminal colors
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

# Run meson introspect and get the key buildoptions as JSON
result = subprocess.run(['meson', 'introspect', 'build', '--buildoptions'],
    stdout=subprocess.PIPE)

# Process the JSON
vars = json.loads(result.stdout)

# Which variables do we want to see?
grab = ['prefix', 'libdir']
ourvars = {}
allgood = True

# Grab the important values
for var in vars:
    if var['name'] in grab:
        ourvars[var['name']] = var['value']

# Show the values for prefix and libdir
print("\nPRESTO's meson build directory currently has:")
print(f"  prefix = {ourvars['prefix']}")
print(f"  libdir = {ourvars['libdir']}\n")

# Now show the user what they should be setting
print(f"Checking for {bcolors.BOLD}PRESTO{bcolors.ENDC} as an environment variable:", end="")
x = os.environ.get('PRESTO')
if x==None:
    print(f"\n  {bcolors.WARNING}WARNING:{bcolors.ENDC} PRESTO environment variable is not set!")
    print("  You can set it using something like:")
    print(f"  export PRESTO={os.getcwd()}")
    allgood = False
else:
    if (os.getcwd()!=os.path.realpath(x)):
        print(f"\n  {bcolors.WARNING}WARNING:{bcolors.ENDC} PRESTO is set, but not to the current directory!")
        allgood = False
    else:
        print(f"  {bcolors.OKGREEN}yes{bcolors.ENDC}")

print(f"Checking for {bcolors.BOLD}TEMPO{bcolors.ENDC} as an environment variable:", end="")
x = os.environ.get('TEMPO')
if x==None:
    print(f"\n  {bcolors.WARNING}WARNING:{bcolors.ENDC} TEMPO environment variable is not set!")
    print("  You need to set it to the location of the TEMPO code.")
    allgood = False
else:
    print(f"  {bcolors.OKGREEN}yes{bcolors.ENDC}")

libinstall = os.path.join(ourvars['prefix'], ourvars['libdir'])
print(f"Is {bcolors.BOLD}{libinstall}{bcolors.ENDC} in LIBRARY_PATH?", end="")
if libinstall not in os.environ.get('LIBRARY_PATH').split(":"):
    print(f"\n  {bcolors.WARNING}WARNING:{bcolors.ENDC} {libinstall} is not in LIBRARY_PATH!")
    print(f"  That may cause issues with linking if {libinstall} is not a standard library directory.")
    print("  You can ensure that it is used by setting LIBRARY_PATH using something like:")
    print(f"  export LIBRARY_PATH={libinstall}:$LIBRARY_PATH")
    allgood = False
else:
    print(f"  {bcolors.OKGREEN}yes{bcolors.ENDC}")

print(f"Is {bcolors.BOLD}{libinstall}{bcolors.ENDC} in LD_LIBRARY_PATH?", end="")
if libinstall not in os.environ.get('LD_LIBRARY_PATH').split(":"):
    print(f"\n  {bcolors.WARNING}WARNING:{bcolors.ENDC} {libinstall} is not in LD_LIBRARY_PATH!")
    print(f"  That may cause issues with running codes if {libinstall} is not a standard library directory.")
    print("  You can ensure that it is used by setting LD_LIBRARY_PATH using something like:")
    print(f"  export LD_LIBRARY_PATH={libinstall}:$LD_LIBRARY_PATH")
    allgood = False
else:
    print(f"  {bcolors.OKGREEN}yes{bcolors.ENDC}")

if allgood:
    print(f"\n{bcolors.OKGREEN}Everything looks good! Let's try to build.{bcolors.ENDC}")
else:
    print(f"\n{bcolors.FAIL}There seem to be some issues. Please fix before trying to build!{bcolors.ENDC}")
    
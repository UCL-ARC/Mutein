import sys
import os
import random
import shutil
import getpass
import pyaes
import base64
import json

def bytes2int(bytes):
    val = 0
    for x in bytes:
        val = (val << 8) + int(x)
    return val

def encrypt_json(args):
    plain_text = getpass.getpass(prompt="Enter JSON as a single line (blank to cancel): ").strip()

    if plain_text == '':
        print("Cancelled")
        return

    initial_counter = os.urandom(16)
    initial_counter_int = bytes2int(initial_counter)
    key_256 = base64.b64decode(os.environ['MUT_PASSWORD'])
    counter = pyaes.Counter(initial_value = initial_counter_int)
    aes = pyaes.AESModeOfOperationCTR(key_256, counter = counter)
    cipher_text = aes.encrypt(plain_text.encode("utf8"))

    message = initial_counter + cipher_text
    with open('./config/password_vault','wb') as f: f.write(message)

def decrypt_json():
    with open('./config/password_vault','rb') as f: message = f.read()

    initial_counter_int = bytes2int(message[:16])
    cipher_text = message[16:]
    key_256 = base64.b64decode(os.environ['MUT_PASSWORD'])
    counter = pyaes.Counter(initial_value = initial_counter_int)
    aes = pyaes.AESModeOfOperationCTR(key_256, counter = counter)
    plain_text = aes.encrypt(cipher_text).decode("utf8")

    return json.loads(plain_text)

def set_password(args):
    password = getpass.getpass(prompt="Enter password (blank for none): ").strip()
    
    if password == '':
        print("")
    else:
        print(password)

def check_cwd(args):
    if os.path.realpath(os.getcwd()) != os.path.realpath(os.environ['MUT_DATA']):
        print(f"Warning: not in expected data directory: {os.environ['MUT_DATA']}")

    exit(0)

def recycle(args):
    recycle_bin = "recycle_bin"

    assert sys.argv[1] == "recycle"
    assert sys.argv[-1] == "END", f"the literal END must follow the file to be recycled: {sys.argv[2:]}"
    assert "END" not in sys.argv[2:-1], f"the literal END must only appear once at the end of the list: {sys.argv[2:]}"

    for fname in sys.argv[2:-1]:
        #sliently ignore non existing files
        if not os.path.exists(fname):
            #not an error if file does not exist
            #this is to prevent scripts having to test for the file before recycling it
            #ie we want any stale file out of the way without using risky "rm -f" command
            continue

        #try to create recycle bin if not already present
        if not os.path.isdir(recycle_bin):
            #fail here if we cannot create the folder
            os.mkdir(recycle_bin)

        new_name = os.path.join(recycle_bin,os.path.basename(fname))

        #make sure not to overwrite existing recycle bin file
        if os.path.exists(new_name):
            new_name += ".%d"%random.randrange(1000000)
            if os.path.exists(new_name): raise Exception(f"File {new_name} already in recycle bin")

        #move to recycle bin
        shutil.move(fname,new_name)

def truncate(args):
    'truncate to zero length without changing the mtime'
    assert sys.argv[1] == "truncate"

    for path in args.file:
        #verify path exists
        if args.n == False: assert os.path.exists(path)

        #get mtime and atime
        mtime = os.path.getmtime(path)
        atime = os.path.getatime(path)

        #truncate file
        f = open(path,'w')
        f.close()

        #reset mtime and atime
        os.utime(path, times=(atime, mtime))

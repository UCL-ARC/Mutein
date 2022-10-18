import sys
import os
import random
import shutil
import getpass
import pyaes
import hashlib
import base64
import json
import csv
#import subprocess
import stat

#salt for password based key derivation funtion
#to prevent reverse lookup of unsalted hash function
#do not change otherwise the passphrase will map to a different key!
salt = base64.b64decode(b'oHP7EUhOZJ37fbTbE/VLbprAEqzSwdHQ2R3FW6/f6xE=')

def encrypt_key(args):
    '''
    prompts for a master password/passphrase interactively from keyboard
    creates a random bulk data password
    encrypts it with the master password
    writes to private file
    along with a random initial AES counter value

    use this command once the setup a key for encfs
    then use decrypt_key at the start of each session to unlock the key
    then delete the plaintext key once done
    '''

    #user should copy-paste in from a password manager
    #or use a long but memorisable passphrase
    master_password = getpass.getpass(prompt="Enter master password: ").strip()
    if master_password == '':
        print("Cancelled")
        return

    #generate random initial counter value for AES counter mode
    #https://github.com/ricmoo/pyaes/blob/master/README.md
    iv = os.urandom(16)

    #derive 256 bit key from the master password
    #https://docs.python.org/3/library/hashlib.html
    master_key = hashlib.pbkdf2_hmac('sha256', master_password.encode('utf8'), salt, 500000, 32)

    #set up AES-256 with IV in counter mode
    aes = pyaes.AESModeOfOperationCTR(master_key,
            counter=pyaes.Counter(initial_value=int.from_bytes(iv, "little")))

    #generate 256 bit bulk data key and encrypt it
    bulk_key = os.urandom(32)

    if args.debug: print('bulk key',base64.b64encode(bulk_key),file=sys.stderr)

    #set to only be readable by user
    f = open(args.output_file,'wb')
    f.close()
    os.chmod(args.output_file,stat.S_IRUSR|stat.S_IWUSR)

    with open(args.output_file,'wb') as f: f.write(iv + aes.encrypt(bulk_key))

def decrypt_key(args):
    '''
    prompts for a master password/passphrase interactively from keyboard
    decrypts a bulk data key from file
    writes to a private file in plaintext
    delete the keyfile once you want to close the session
    set up the initial encrypted key using encrypt_key
    '''

    #user should copy-paste in from a password manager
    #or use a long but memorisable passphrase
    master_password = getpass.getpass(prompt="Enter master password: ").strip()
    if master_password == '':
        print("Cancelled")
        return

    #derive 256 bit key from the master password
    #https://docs.python.org/3/library/hashlib.html
    master_key = hashlib.pbkdf2_hmac('sha256', master_password.encode('utf8'), salt, 500000, 32)

    with open(args.input_file,'rb') as f: data = f.read()
    iv = data[:16]
    encrypted_key = data[16:]

    #set up AES-256 with IV in counter mode
    aes = pyaes.AESModeOfOperationCTR(master_key,
            counter=pyaes.Counter(initial_value=int.from_bytes(iv, "little")))

    #decrypt 256 bit bulk data key to file
    bulk_key = aes.decrypt(encrypted_key)

    if args.debug: print('bulk key',base64.b64encode(bulk_key),file=sys.stderr)

    #set to only be readable by user
    f = open(args.output_file,'wb')
    f.close()
    os.chmod(args.output_file,stat.S_IRUSR|stat.S_IWUSR)

    with open(args.output_file,'wb') as f: f.write(bulk_key)

# def get_key(args):
#     if args.key_file == "PROMPT":
#         key = getpass.getpass(prompt="Enter password (blank to cancel): ").strip()
#         if key == '':
#             print("Cancelled")
#             exit(0)
#     else:
#         key = open(args.key_file).read().strip()
#     return key

# def aes_encrypt(args):
#     '''
#     AES encrypt of a bulk a file
#     '''
#     key = get_key(args)
#     cmd = "gpg --symmetric --batch "
#     cmd += f"--passphrase {key} --cipher-algo AES128 --output {args.output_file} {args.input_file}"
#     subprocess.run(cmd,shell=False,check=True)

# def aes_decrypt(args):
#     '''
#     AES decrypt of a bulk file
#     '''
#     key = get_key(args)
#     cmd = "gpg --decrypt --batch "
#     cmd += f"--passphrase {key} --cipher-algo AES128 --output {args.output_file} {args.input_file}"
#     subprocess.run(cmd,shell=False,check=True)

def symlink(args):
    '''
    create a relative symlink with the shortest path to the target
    regardless of the current working directory
    better than doing eg:
    ln -sr datasets/set1/subset1/accession1/file datasets/set1/subset1/accession1/link
    '''
    target_dir = os.path.dirname(args.target)
    target_file = os.path.basename(args.target)
    link_dir = os.path.dirname(args.link)

    #find path to target relative to symlink
    relpath = os.path.relpath(target_dir,start=link_dir)

    final_path = os.path.join(relpath,target_file)

    if os.path.exists(args.link):
        if os.path.islink(args.link):
            #overwrite any existing link
            os.unlink(args.link)
        else:
            #but refuse to overwrite any existing file
            raise Exception(f"path exists but is not a symlink: {args.link}")

    os.symlink(final_path,args.link)

def cut(args):
    '''
    parse stdin through csv reader using pass-through options
    output to stdout the fields requested from the -f option
    '''

    if args.csvopts != None:
        kwargs = json.loads(args.csvopts)
    else:
        kwargs = {}

    reader = csv.reader(sys.stdin,**kwargs)

    #output field delimiter
    if args.output_delim == None:
        delim = reader.dialect.delimiter
    else:
        delim = args.output_delim

    for row in reader:
        print(delim.join([row[i-1] for i in args.field]))


def bytes2int(bytes):
    val = 0
    for x in bytes:
        val = (val << 8) + int(x)
    return val

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

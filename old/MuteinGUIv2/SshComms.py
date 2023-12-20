from threading import Thread

class AsyncSsh(Thread):
    def __init__(self):
        super().__init__()

        self.user = ""
        self.pw = ""
        self.server = ""
        self.command = ""
        self.msg = ""
        self.error = ""


    def setDetails(self,user,pw,server,command):
        self.user = user
        self.pw = pw
        self.server = server
        self.command = command
        
    def getSshShell(self):
        
        import spur
        try:
            shell = spur.SshShell(hostname= self.server + ".rc.ucl.ac.uk", username=self.user, password=self.pw)
            return shell
        except:
            return None
    
    
    def run(self):
        try:
            shell = self.getSshShell()
            result = shell.run(self.command)
            self.msg = result.output  
            self.error = result.stderr_output      
        except:
            self.error = "Some kind of error"

   
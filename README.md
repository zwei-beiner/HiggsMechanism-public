## How to build from source (On a Mac)
1. **Clone the git repository**  
```$ git clone https://github.com/zwei-beiner/HiggsMechanism.git```  
This creates a new directory ```./HiggsMechanism``` .  
  


2. **Navigate into the new directory**  
```$ cd ./HiggsMechanism```  


3. **Create a new python virtual environment**  
```$ python3 -m venv venv```


4. **Activate the virtual environment**  
```$ source venv/bin/activate```


5. **Update pip**  
```(venv)$ pip install --upgrade pip setuptools```  


6. **Install the required dependencies**  
```(venv)$ pip install PyQt5 numpy pyqtgraph pyinstaller```  


7. **Run the build command**  
```(venv)$ pyinstaller --name "Higgs Simulator" --icon icon.ico --windowed main.py```  
This creates the directories ```build``` and  ```dist``` and the file ```Higgs Simulator.spec```.  
The .app executable can be found in ```dist```. Additionally, a command-line executable can be found at ```dist/main/main```.

## If building from source goes wrong
Standard practice on using pyinstaller documentation and StackOverflow applies.  

To get more debug information when building, use the ```--debug=all``` flag:  
```(venv)$ pyinstaller --debug=all main.py```  

(Note: If ```pyinstaller``` gives a ```WARNING``` about missing imports, they can possibly be ignored, e.g. warnings about ```pkg_resources.py2_warn``` and ```pkg_resources.markers``` can be ignored at the time of writing. To test which missing imports can be ignored, write a small demo program which is guaranteed to work and run pyinstaller on it, then look at which missing imports pyinstaller is complaining about.)


## How the image icon.ico was made
.ico file was created by uploading a square high-resolution png to 
https://icoconvert.com
and selecting "Custom sizes", "Multi-size in one icon" and selecting all 9 sizes, then click "Convert ICO"
(We have to do this because an .ico file stores all possible resolutions we might want to display, up to 256x256.)


## How to build from source (on Windows)
Microsoft releases free Windows virtual machines for development (with the catch that the Windows license has an expiry date).
These can be downloaded here:  
https://developer.microsoft.com/en-us/windows/downloads/virtual-machines/

VirtualBox was used and the resulting virtual machine was 50GB (!) in size.

In Windows, to build the ".exe" file,

1. **Install git and python.**  
This can be done by downloading the installers from the websites. 
 

2. Follow the instructions above (for "On a Mac") with the changes:
   1. Run all commands in the command line as an administrator. To do this,
      1. Click on "Start" 
      2. Type "cmd"
      3. Right-click the icon that shows up and click "As an administrator"
   2. "python3" becomes just "python"
   3. Activating the venv is done with the command ```.\venv\Scripts\activate```
   4. The build command is  

   ```pyinstaller --name "Higgs Simulator" --icon icon.ico --windowed --onefile main.py```

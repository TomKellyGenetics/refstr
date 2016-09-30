##Install Python language
```shell
sudo apt-get install python
```

##Install Python packages
```shell
sudo apt-get install python-pip
sudo pip install numpy
sudo pip install pandas
```

##If headers missing
```shell
sudo apt-get install python-dev
```

For recent linux distros may use:
```shell
sudo apt-get install python3 python3-dev
```

To resolve version conflicts
```shell
sudo apt-get install aptitude
sudo aptitude install python-dev
#no to first option (don't install any package dependencies)
#yes to second option to downgrade packages
```


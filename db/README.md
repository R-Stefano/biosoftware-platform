# Ion-Channel Database.
## A python script to build Ion-Channle Database for Web Platform
The code in this directory is used to build a sqlite file to use in the [EDEn](https://www.sciencedirect.com/science/article/pii/S2589004218302323#:~:text=The%20electroceutical%20design%20environment%20(EDEn,a%20range%20of%20disease%20states)) Web Platform

# Setup
Install virtualenv
```
pip install virtualenvwrapper
mkvirtualenv biosoftware
```

Make executable the bash file
```
chmod +x setup_icdb.sh
```

Activate virtualenv and install libraries
```
workon biosoftware
pip install -r requirements.txt
```

Run Setup
```
./setup_icdb.sh 
```

# Contact
pwinter@ualberta.ca

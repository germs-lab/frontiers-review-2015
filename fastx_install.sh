sudo apt-get -y install gcc g++ pkg-config
cd libgtextutils-0.7
sudo ./configure
sudo make
sudo make install
cd ..
cd fastx_toolkit-0.0.14
sudo ./configure
sudo make 
sudo make install
cd ..

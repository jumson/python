# python_show
showing some example python scripts I made for various reasons

* **raw_to_symbol.py** is a script that consumes a raw iq sample file produced by a Software Defined Radio; determines the signal, and outputs symboles (encoded bits) along with measurement information about the signal itself.
* **basic_help.py** is a file required by the script above; I use this like a C/C++ *#include* file --- a library of classes and functions that I used in many other scripts. 
* **csv2bin.py** For Saleae generated analysis (I2C, Async Serial, SPI, etc) CSV files, pulls out the bytes and puts them into a binary file.

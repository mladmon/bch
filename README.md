# BCH error-correcting code implementation in MATLAB

This project was part of a directed study in the EECS department of the University of Michigan in the Winter, 2011 semester. The group I worked with conducted research on memristor memory technologies. I apologize for the lack of documentation; I wrote this code as a junior in undergrad prior to appreciating the importance of proper documentation. However, upon inspecting the code eight years later, I'm proud to say that it seems to be self-documented and easy to follow (assuming you're familiar with BCH theory (which is a very large assumption to make)). The header of each file describes its function, and there are comments throughout the main script, [BCH_Code.m](https://github.com/mladmon/bch/blob/master/BCH_Code.m).

This MATLAB code implements a (31, 16) BCH code that's capable of correcting up to 3 errors:
   - 31-bit codeword
   - 16-bit message
   - 3 correctable errors

## Running the script
To run the script, simply run the following command in the MATLAB interpreter:

```
>> BCH_Code
```

You should get the following output:

![alt text](bch_code.png)

This example demonstrates how a 31-bit codeword endures 3 errors during transmission, and how these errors are subsequently detected and corrected in the received codeword. Note, you can change the bits that get flipped during transmission by modifying the assignments on lines 46-48. For example, try changing them to:

```
rx_codeword(1) = 1;
rx_codeword(3) = 1;
rx_codeword(7) = 1;
```

and run the script. Have fun!

# No-Reference PSNR Estimation Algorithm for H.264 Encoded Video Sequences 

This repository contains a no-reference PSNR calculation algorithm
implementation in C.

The program takes in input a single XML file and outputs several info in
different files. 

The input file has to be in XML form, see the attached `partialTmpDebug.xml`
for reference. In order to create such a file from a video source, a proper
decoder has to be used. 
For a `.h264` source file, the `ldecod_Ghent_XML.exe` program has been used.
You may find it [here](http://vqegstl.ugent.be).
This C implementation extracts from the XML input file just the information 
needed for the calculations. 

The output files are the following:
* <output_file_name>.frameMseBeta
* <output_file_name>.frameMseLambda
* <output_file_name>.framePsnrBeta
* <output_file_name>.framePsnrLambda
* <output_file_name>.mseBeta
* <output_file_name>.mseLambda
* <output_file_name>.psnrBeta
* <output_file_name>.psnrLambda

Where the first 4 have one value for each Picture whilst the last 4 have
this structure:
```
mse_MB1,mse2_MB,...,mse_MBN\n
mse_MB1,mse_MB2,...,mse_MBN\n
...
mse_MB1,mse_MB2,...,mse_MBN
```
(same for PSNR).

# Launch
Launch by using the internal makefile. The makefile contains the following:
```
$ gcc no_ref_psnr.c -Wall -lm -o no_ref_psnr 
$ ./no_ref_psnr input_file_name.xml output_file_name 
```

# Debug 
This repository contains a `.c` file and a `.h` file attached.
For debug purposes, a dedicated `DEBUG` flag can be set to `1` in the header
file.
Remember to compile with the `-g` flag to debug with gdb or valgrind. 

# Credits
This implementation is based on the paper entitled:
"No-ref psnr estimation algorithm"
by Brandao and Queluz. 

You are free to use this program (see LICENSE) but if you do so please include
the paper called:
```
M. Siekkinen, T. Kamarainen, L. Favario, E. Masala, "Can You See What I See?
Quality-of-Experience Measurements of Mobile Live Video Broadcasting", ACM
Transactions on Multimedia Computing, Communications, and Applications (TOMM),
Volume 14 Issue 2s, May 2018 
```
in the references section of your papers.
You may find it [here](https://dx.doi.org/10.1145/3165279)

# License
This program is covered by a GNU GPLv3 license. See the `LICENSE.md` file for
more extended details. 

## Disclaimer
```
This software is provided by the copyright holders and contributors “AS IS” and
any express or implied warranties, including, but not limited to, the implied
WARRANTIES of MERCHANTABILITY and FITNESS for a particular purpose are
DISCLAIMED. In NO EVENT shall the copyright owner or contributors be liable for
any DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, or CONSEQUENTIAL DAMAGES
(including, but not limited to, procurement of substitute goods or services;
loss of use, data, or profits; or business interruption) however caused and on
any theory of LIABILITY, whether in contract, strict liability, or tort
(including negligence or otherwise) arising in ANY WAY out of the use of this
software, even if advised of the possibility of such damage.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
```

# RNS-AVX-POC

A Proof-of-Concept bignum library for cryptography using Remainder Number Systems and Montgomery modular multiplitcation algorithm. Utilizes AVX2 for acceleration.

This PoC rovides a 2048-bit unsigned int type and a set of arithmetic operations that could be used to implement some public-key cryptographic algorithms. An example implementation of RSA is provided in rsa.c. Types of other lengths could be implemented using the same algorithms. 

Acknowledgments:

1. E. Ochoa-Jiménez, L. Rivera-Zamarripa, N. Cruz-Cortés and F. Rodríguez-Henríquez, "Implementation of RSA Signatures on GPU and CPU Architectures," in IEEE Access, vol. 8, pp. 9928-9941, 2020, doi: 10.1109/ACCESS.2019.2963826.
   
2. Kawamura, S., Koike, M., Sano, F., Shimbo, A. (2000). Cox-Rower Architecture for Fast Parallel Montgomery Multiplication. In: Preneel, B. (eds) Advances in Cryptology — EUROCRYPT 2000. EUROCRYPT 2000. Lecture Notes in Computer Science, vol 1807. Springer, Berlin, Heidelberg. https://doi.org/10.1007/3-540-45539-6_37

3. S. Kawamura, Y. Komano, H. Shimizu, T. Yonemura. “RNS Montgomery reduction algorithms using quadratic residuosity”, Journal of Cryptographic Engineering, 9:4 (2019) 

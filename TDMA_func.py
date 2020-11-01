# -*- coding: utf-8 -*-

'''
Cette fonction effectue la resolution d'un systeme d'equation lineaire
[T][X] = [C]
Methode presente dans le livre de Anderson,1995, Appendix A: Thomas algorithm for the solution of a tridiagonal system of equation
ou la matrice des coefficients, [T] est une matrice de format n x n Tridiagonale
ou les elements de la diagonale D sont plus grands que tous les autres
elements sur la meme ligne (c'est a dire que l'element a(1,1) est le plus
grand elements de la ligne 1, l'element a(2,2) le plus grand element de la
ligne 2 et ainsi de suite.
--------------------------------------------------------------------------
Les variables d'entrees sonts
--------------------------------------------------------------------------
Les elements de la diagonale inferieure sous forme de vecteur :
B = [a(2,1);a(3,2);a(4,3);...;a(n,n-1)]
Les elements de la diagonale sous forme de vecteur:
D = [a(1,1);a(2,2);a(3,3);...;a(n,n)]
Les elements de la diagonale superieure sous forme de vecteur:
A = [a(1,2);a(2,3);a(3,4);...;a(n,n+1)]
Les elements du vecteur de constante :
C = [b(1,1);b(2,1);b(3,1);...;b(n,1)]

-------------------------------------------------------------------------
Les variables de sortie sont
--------------------------------------------------------------------------
X le vecteur solution
'''
# Library
import numpy as np

def TDMA(A,B,C,D):
    
    #==========================================================================
    n = D.size-1; # nombre d'inconues (ou de lignes)
    # Attention a la longueur des vecteurs A et B si A, B, C, D sont de la meme longueur commente les ligne 30 et 31!
    B = np.insert(B,0,0);
    A = np.append(A,0);
    X = np.zeros([n+1,1]); # initialisation du vecteur solution
    
    # Obtention d'une matrice bidiagonale superieure
    for i in np.arange(0,n+1,1):
        D[i] = D[i]-(B[i]*A[i-1])/D[i-1]; # On remplace les coefficents sur la diagonale principale
        C[i] = C[i]-(C[i-1]*B[i])/D[i-1]; # On remplace les termes dans le membre de droite 
    
    # resolution du dernier element du vecteur solution
    X[n] = C[n]/D[n];
    
    # Resolution de tous les autres elements du vecteur solution en travaillant
    # a reculons vers le haut
    
    for i in np.arange(n-1,-1,-1):
        X[i] = (C[i]-A[i]*X[i+1])/D[i];
        
    return X


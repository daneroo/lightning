// Decompiled by Jad v1.5.8e. Copyright 2001 Pavel Kouznetsov.
// Jad home page: http://www.geocities.com/kpdus/jad.html
// Decompiler options: packimports(3) 
// Source File Name:   Complex.java

import java.io.PrintStream;

public class Complex
{

    public Complex(double d, double d1)
    {
        re = d;
        im = d1;
    }

    public String toString()
    {
        return re + " + " + im + "i";
    }

    public double abs()
    {
        return Math.sqrt(re * re + im * im);
    }

    public Complex plus(Complex complex)
    {
        Complex complex1 = this;
        double d = complex1.re + complex.re;
        double d1 = complex1.im + complex.im;
        Complex complex2 = new Complex(d, d1);
        return complex2;
    }

    public Complex minus(Complex complex)
    {
        Complex complex1 = this;
        double d = complex1.re - complex.re;
        double d1 = complex1.im - complex.im;
        Complex complex2 = new Complex(d, d1);
        return complex2;
    }

    public Complex times(Complex complex)
    {
        Complex complex1 = this;
        double d = complex1.re * complex.re - complex1.im * complex.im;
        double d1 = complex1.re * complex.im + complex1.im * complex.re;
        Complex complex2 = new Complex(d, d1);
        return complex2;
    }

    public Complex times(double d)
    {
        return new Complex(d * re, d * im);
    }

    public Complex conjugate()
    {
        return new Complex(re, -im);
    }

    public static Complex plus(Complex complex, Complex complex1)
    {
        double d = complex.re + complex1.re;
        double d1 = complex.im + complex1.im;
        Complex complex2 = new Complex(d, d1);
        return complex2;
    }

    public static void main(String args[])
    {
        Complex complex = new Complex(5D, 6D);
        System.out.println("a = " + complex);
        Complex complex1 = new Complex(-3D, 4D);
        System.out.println("b = " + complex1);
        Complex complex2 = complex1.times(complex);
        System.out.println("c = " + complex2);
        Complex complex3 = complex2.conjugate();
        System.out.println("d = " + complex3);
        double d = complex1.abs();
        System.out.println("e = " + d);
        Complex complex4 = complex.plus(complex1);
        System.out.println("f = " + complex4);
        Complex complex5 = plus(complex, complex1);
        System.out.println("g = " + complex5);
    }

    public final double re;
    public final double im;
}

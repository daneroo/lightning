// Decompiled by Jad v1.5.8e. Copyright 2001 Pavel Kouznetsov.
// Jad home page: http://www.geocities.com/kpdus/jad.html
// Decompiler options: packimports(3) 
// Source File Name:   FFT.java

import java.awt.Graphics2D;
import java.awt.image.BufferedImage;

public class FFT
{

    public FFT(BufferedImage bufferedimage)
    {
        filterWidth = bufferedimage.getWidth();
        filterHeight = bufferedimage.getHeight();
        int ai[] = bufferedimage.getRGB(0, 0, filterWidth, filterHeight, null, 0, filterWidth);
        filter = new Complex[filterWidth][];
        for(int i = 0; i < filterHeight; i++)
        {
            filter[i] = new Complex[filterHeight];
            for(int j = 0; j < filterWidth; j++)
            {
                int k = j + i * filterWidth;
                int l = ai[k];
                int i1 = l & 0xff;
                filter[i][j] = new Complex((float)i1 / 255F, 0.0D);
            }

        }

        fft2D(filter, filterWidth, filterHeight);
        rendering = null;
    }

    public void fft2D(Complex acomplex[][], int i, int j)
    {
        for(int k = 0; k < j; k++)
            acomplex[k] = fft(acomplex[k]);

        for(int l = 0; l < i; l++)
        {
            Complex acomplex1[] = new Complex[j];
            for(int i1 = 0; i1 < j; i1++)
                acomplex1[i1] = acomplex[i1][l];

            Complex acomplex2[] = fft(acomplex1);
            for(int j1 = 0; j1 < j; j1++)
                acomplex[j1][l] = acomplex2[j1];

        }

    }

    public void ifft2D(Complex acomplex[][], int i, int j)
    {
        for(int k = 0; k < j; k++)
            acomplex[k] = ifft(acomplex[k]);

        for(int l = 0; l < i; l++)
        {
            Complex acomplex1[] = new Complex[j];
            for(int i1 = 0; i1 < j; i1++)
                acomplex1[i1] = acomplex[i1][l];

            Complex acomplex2[] = ifft(acomplex1);
            for(int j1 = 0; j1 < j; j1++)
                acomplex[j1][l] = acomplex2[j1];

        }

    }

    public BufferedImage convolve2D(BufferedImage bufferedimage, int ai[])
    {
        rendering = new float[filterWidth * filterHeight];
        int ai1[] = bufferedimage.getRGB(0, 0, bufferedimage.getWidth(), bufferedimage.getHeight(), null, 0, bufferedimage.getWidth());
        int i = bufferedimage.getWidth();
        int j = bufferedimage.getHeight();
        image = new Complex[i][];
        for(int j1 = 0; j1 < j; j1++)
        {
            image[j1] = new Complex[j];
            for(int k = 0; k < i; k++)
            {
                int i2 = k + j1 * i;
                int j2 = ai1[i2];
                int k2 = j2 & 0xff;
                image[j1][k] = new Complex((float)k2 / 255F, 0.0D);
            }

        }

        fft2D(image, i, j);
        for(int k1 = 0; k1 < j; k1++)
        {
            for(int l = 0; l < i; l++)
                image[k1][l] = image[k1][l].times(filter[k1][l]);

        }

        ifft2D(image, i, j);
        for(int l1 = 0; l1 < j; l1++)
        {
            for(int i1 = 0; i1 < i; i1++)
            {
                float f = (float)image[l1][i1].re;
                rendering[i1 + l1 * i] = f;
            }

        }

        return relight(ai);
    }

    public BufferedImage relight(int ai[])
    {
        if(rendering == null)
            return null;
        BufferedImage bufferedimage = new BufferedImage(filterWidth, filterHeight, 1);
        Graphics2D graphics2d = (Graphics2D)bufferedimage.getGraphics();
        for(int i = 0; i < filterHeight; i++)
        {
            for(int j = 0; j < filterWidth; j++)
            {
                float f = rendering[j + i * filterWidth];
                float f1 = (f * (float)ai[0]) / 255F;
                float f2 = (f * (float)ai[1]) / 255F;
                float f3 = (f * (float)ai[2]) / 255F;
                f1 = f1 <= 1.0F ? f1 : 1.0F;
                f2 = f2 <= 1.0F ? f2 : 1.0F;
                f3 = f3 <= 1.0F ? f3 : 1.0F;
                int k = ((int)(f1 * 255F) << 16) + ((int)(f2 * 255F) << 8) + (int)(f3 * 255F);
                bufferedimage.setRGB(j, i, k);
            }

        }

        return bufferedimage;
    }

    public static Complex[] fft(Complex acomplex[])
    {
        int i = acomplex.length;
        Complex acomplex1[] = new Complex[i];
        if(i == 1)
        {
            acomplex1[0] = acomplex[0];
            return acomplex1;
        }
        if(i % 2 != 0)
            throw new RuntimeException("N is not a power of 2");
        Complex acomplex2[] = new Complex[i / 2];
        Complex acomplex3[] = new Complex[i / 2];
        for(int j = 0; j < i / 2; j++)
            acomplex2[j] = acomplex[2 * j];

        for(int k = 0; k < i / 2; k++)
            acomplex3[k] = acomplex[2 * k + 1];

        Complex acomplex4[] = fft(acomplex2);
        Complex acomplex5[] = fft(acomplex3);
        for(int l = 0; l < i / 2; l++)
        {
            double d = ((double)(-2 * l) * 3.1415926535897931D) / (double)i;
            Complex complex = new Complex(Math.cos(d), Math.sin(d));
            acomplex1[l] = acomplex4[l].plus(complex.times(acomplex5[l]));
            acomplex1[l + i / 2] = acomplex4[l].minus(complex.times(acomplex5[l]));
        }

        return acomplex1;
    }

    public static Complex[] ifft(Complex acomplex[])
    {
        int i = acomplex.length;
        for(int j = 0; j < i; j++)
            acomplex[j] = acomplex[j].conjugate();

        Complex acomplex1[] = fft(acomplex);
        for(int k = 0; k < i; k++)
            acomplex1[k] = acomplex1[k].conjugate();

        for(int l = 0; l < i; l++)
            acomplex1[l] = acomplex1[l].times(1.0D / (double)i);

        return acomplex1;
    }

    public static Complex[] convolve(Complex acomplex[], Complex acomplex1[])
    {
        if(acomplex.length != acomplex1.length)
            throw new RuntimeException("Dimensions don't agree");
        int i = acomplex.length;
        Complex acomplex2[] = fft(acomplex);
        Complex acomplex3[] = fft(acomplex1);
        Complex acomplex4[] = new Complex[i];
        for(int j = 0; j < i; j++)
            acomplex4[j] = acomplex2[j].times(acomplex3[j]);

        return ifft(acomplex4);
    }

    public int filterWidth;
    public int filterHeight;
    Complex filter[][];
    Complex image[][];
    float rendering[];
}

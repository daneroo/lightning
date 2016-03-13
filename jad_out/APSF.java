// Decompiled by Jad v1.5.8e. Copyright 2001 Pavel Kouznetsov.
// Jad home page: http://www.geocities.com/kpdus/jad.html
// Decompiler options: packimports(3) 
// Source File Name:   APSF.java

import java.util.Vector;

public class APSF
{

    public float[] getKernel()
    {
        return _kernel;
    }

    public APSF(int i)
    {
        _res = i;
        if(_res % 2 == 0)
            _res++;
        _kernel = new float[_res * _res];
        _kernel1D = new float[_res];
        _q = 0.99F;
        _R = 200F;
        _D = 2000F;
        _T = 1.001F;
        _sigma = 0.5F;
        _maxTerms = 200;
        _I0 = 1.0F;
        _retinaSize = 0.01F;
        _eyeSize = 0.025F;
    }

    private float gM(float f, int i)
    {
        return i != 0 ? (float)Math.exp(-((double)(betaM(i, _q) * _T) + (double)alphaM(i) * Math.log(_T))) : 0.0F;
    }

    private float alphaM(float f)
    {
        return f + 1.0F;
    }

    private float betaM(float f, float f1)
    {
        return ((2.0F * f + 1.0F) / f) * (1.0F - (float)Math.pow(f1, (int)f - 1));
    }

    private float factorial(float f)
    {
        return f > 1.0F ? f * factorial(f - 1.0F) : 1.0F;
    }

    private float choose(float f, float f1)
    {
        return factorial(f) / (factorial(f1) * factorial(f - f1));
    }

    private float legendreM(int i, float f)
    {
        Vector vector = new Vector();
        vector.addElement(new Float(1.0F));
        vector.addElement(new Float(f));
        for(int j = 2; j <= i; j++)
        {
            float f1 = ((Float)vector.elementAt(j - 1)).floatValue();
            float f2 = ((Float)vector.elementAt(j - 2)).floatValue();
            float f3 = ((2.0F * (float)j - 1.0F) * f * f1 - ((float)j - 1.0F) * f2) / (float)j;
            vector.addElement(new Float(f3));
        }

        return ((Float)vector.elementAt(i)).floatValue();
    }

    private float pointAPSF(float f)
    {
        float f1 = 0.0F;
        for(int i = 0; i < _maxTerms; i++)
            f1 += (gM(_I0, i) + gM(_I0, i + 1)) * legendreM(i, f);

        return f1;
    }

    public void generateKernel()
    {
        float f = _retinaSize / (float)_res;
        float f1 = _retinaSize / (float)_res;
        int i = _res / 2;
        int j = 0;
        float f2 = 0.0F;
        float f3 = 1000F;
        for(int k = 0; k < _res; k++)
        {
            for(int i1 = 0; i1 < _res;)
            {
                float f5 = (float)(i1 - i) * f;
                float f6 = (float)(k - i) * f1;
                float f7 = (float)Math.sqrt(f5 * f5 + f6 * f6);
                if(f7 / _eyeSize > _R / _D)
                {
                    _kernel[j] = 0.0F;
                } else
                {
                    float f8 = -f7 * f7 * _D * _D + _eyeSize * _eyeSize * _R * _R + f7 * f7 * _R * _R;
                    f8 = _eyeSize * _eyeSize * _D - _eyeSize * (float)Math.sqrt(f8);
                    f8 /= _eyeSize * _eyeSize + f7 * f7;
                    float f9 = 3.141593F - (float)Math.atan(_retinaSize / f7) - (float)Math.asin((_D - f8) / _R);
                    _kernel[j] = pointAPSF((float)Math.cos(f9));
                    f3 = _kernel[j] >= f3 ? f3 : _kernel[j];
                }
                f2 = _kernel[j] <= f2 ? f2 : _kernel[j];
                i1++;
                j++;
            }

        }

        if(f3 > 0.0F)
        {
            for(int l = 0; l < _res * _res; l++)
                if(_kernel[l] > 0.0F)
                    _kernel[l] -= f3;

            f2 -= f3;
        }
        if(f2 > 1.0F)
        {
            float f4 = 1.0F / f2;
            for(int j1 = 0; j1 < _res * _res; j1++)
                _kernel[j1] *= f4;

        }
    }

    public void generateKernel1D()
    {
        float f = _retinaSize / (float)_res;
        float f1 = _retinaSize / (float)_res;
        int i = _res / 2;
        float f2 = 0.0F;
        float f3 = 1000F;
        float f4 = 0.0F;
        for(int j = 0; j < _res; j++)
        {
            float f6 = (float)(j - i) * f;
            float f7 = (float)Math.sqrt(f6 * f6 + f4 * f4);
            if(f7 / _eyeSize > _R / _D)
            {
                _kernel1D[j] = 0.0F;
            } else
            {
                float f8 = -f7 * f7 * _D * _D + _eyeSize * _eyeSize * _R * _R + f7 * f7 * _R * _R;
                f8 = _eyeSize * _eyeSize * _D - _eyeSize * (float)Math.sqrt(f8);
                f8 /= _eyeSize * _eyeSize + f7 * f7;
                float f9 = 3.141593F - (float)Math.atan(_retinaSize / f7) - (float)Math.asin((_D - f8) / _R);
                _kernel1D[j] = pointAPSF((float)Math.cos(f9));
                f3 = _kernel1D[j] >= f3 ? f3 : _kernel1D[j];
            }
            f2 = _kernel1D[j] <= f2 ? f2 : _kernel1D[j];
        }

        if(f3 > 0.0F)
        {
            for(int k = 0; k < _res; k++)
                if(_kernel1D[k] > 0.0F)
                    _kernel1D[k] -= f3;

            f2 -= f3;
        }
        if(f2 > 1.0F)
        {
            float f5 = 1.0F / f2;
            for(int l = 0; l < _res; l++)
                _kernel1D[l] *= f5;

        }
    }

    float _q;
    float _T;
    float _I0;
    float _sigma;
    float _R;
    float _D;
    float _retinaSize;
    float _eyeSize;
    int _maxTerms;
    int _res;
    public float _kernel[];
    public float _kernel1D[];
}

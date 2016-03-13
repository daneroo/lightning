// Decompiled by Jad v1.5.8e. Copyright 2001 Pavel Kouznetsov.
// Jad home page: http://www.geocities.com/kpdus/jad.html
// Decompiler options: packimports(3) 
// Source File Name:   Lightning.java

import java.awt.*;
import java.awt.geom.Line2D;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.util.*;

public class Lightning
    implements Runnable
{

    public int xRes()
    {
        return _xRes;
    }

    public int yRes()
    {
        return _yRes;
    }

    public float[] phase()
    {
        return _phase;
    }

    public float[] potential()
    {
        return _potential;
    }

    public float eta()
    {
        return _power;
    }

    public void setEta(float f)
    {
        _power = f;
    }

    public Lightning(int i, int j, float f, BackBuffer backbuffer, TextArea textarea)
    {
        backBuffer = backbuffer;
        statusBox = textarea;
        Date date = new Date();
        long l = date.getHours() * 3600 + date.getMinutes() * 60 + date.getSeconds();
        rand = new Random(l);
        init(i, j, f);
        pauseFlag = true;
        apsf = null;
        longestBranch = 0.05F;
        chargeType = -1;
    }

    public void init(int i, int j, float f)
    {
        _xRes = i;
        _yRes = j;
        _totalCells = _xRes * _yRes;
        _width = 0.015625F * (float)i;
        _height = 0.015625F * (float)j;
        _dx = _width / (float)_xRes;
        _dy = _height / (float)_yRes;
        _phase = new float[_totalCells];
        _potential = new float[_totalCells];
        _charge = new float[_totalCells];
        _lightning = new boolean[_totalCells];
        _residual = new float[_totalCells];
        _direction = new float[_totalCells];
        _q = new float[_totalCells];
        _s = new float[_totalCells];
        root = null;
        DAGHash = new HashMap();
        nodeOrder = new Vector();
        _power = 1.0F;
        initPotential();
        leaderCharge = 0.75F;
        secondaryCharge = 0.3F;
        sideCharge = 0.1F;
        drawGridAndTree();
        Date date = new Date();
        long l = date.getHours() * 3600 + date.getMinutes() * 60 + date.getSeconds();
        rand = new Random(l);
    }

    public void initPotential()
    {
        for(int i = 0; i < _totalCells; i++)
        {
            _phase[i] = _potential[i] = _residual[i] = _direction[i] = _q[i] = _s[i] = _charge[i] = 0.0F;
            _lightning[i] = false;
        }

        for(int j = 0; j < _xRes; j++)
            _potential[j] = 1.0F;

        for(int l = (int)((float)_yRes * 0.9F); l < _yRes; l++)
        {
            for(int k = 0; k < _xRes; k++)
            {
                if(k != _xRes / 2)
                    continue;
                int i1 = k + l * _xRes;
                _potential[i1] = 0.0F;
                _phase[i1] = 1.0F;
                _lightning[i1] = true;
                if(root == null)
                {
                    root = new DAGNode(i1);
                    DAGHash.put(new Integer(i1), root);
                } else
                {
                    addSegment(i1);
                }
            }

        }

    }

    private void addSegment(int i)
    {
        int j = i;
        int k = i / _xRes;
        int l = i - k * _xRes;
        if(l == 0 || k == 0 || l == _xRes - 1 || k == _yRes - 1)
            return;
        if(_lightning[(i - 1) + _xRes])
            j = (i - 1) + _xRes;
        if(_lightning[i + _xRes])
            j = i + _xRes;
        if(_lightning[i + 1 + _xRes])
            j = i + 1 + _xRes;
        if(_lightning[i + 1])
            j = i + 1;
        if(_lightning[(i + 1) - _xRes])
            j = (i + 1) - _xRes;
        if(_lightning[i - _xRes])
            j = i - _xRes;
        if(_lightning[i - 1 - _xRes])
            j = i - 1 - _xRes;
        if(_lightning[i - 1])
            j = i - 1;
        if(j == i)
        {
            return;
        } else
        {
            DAGNode dagnode = (DAGNode)DAGHash.get(new Integer(j));
            DAGNode dagnode1 = new DAGNode(i);
            dagnode1.parent = dagnode;
            dagnode.children.addElement(dagnode1);
            dagnode1.potential = _potential[i];
            DAGHash.put(new Integer(i), dagnode1);
            nodeOrder.addElement(new Integer(i));
            return;
        }
    }

    public void drawGridAndTree()
    {
        BackBuffer _tmp = backBuffer;
        int i = 512;
        BackBuffer _tmp1 = backBuffer;
        int j = 512;
        BufferedImage bufferedimage = new BufferedImage(i, j, 1);
        Graphics2D graphics2d = (Graphics2D)bufferedimage.getGraphics();
        float f = (float)i / (float)_xRes;
        float f1 = (float)j / (float)_yRes;
        for(int k = 0; k < _yRes; k++)
        {
            for(int l = 0; l < _xRes; l++)
            {
                int i1 = l + _xRes * k;
                if(_charge[i1] >= 1.0F)
                    graphics2d.setPaint(new Color(0, 255, 0));
                else
                if(_lightning[i1] && _potential[i1] < 1.0F)
                {
                    graphics2d.setPaint(new Color(0, 0, 255));
                } else
                {
                    int j1 = (int)(_potential[i1] * 255F);
                    graphics2d.setPaint(new Color(j1, 0, 0));
                }
                graphics2d.fill(new java.awt.geom.Rectangle2D.Double(f * (float)l, f1 * (float)(_yRes - 1 - k), f, f1));
            }

        }

        drawSegments(root, graphics2d, i, j, 0, 0);
        backBuffer.backBuffer = bufferedimage;
        backBuffer.backBufferContext = graphics2d;
    }

    public void drawSegments(DAGNode dagnode, Graphics2D graphics2d, int i, int j, int k, int l)
    {
        int i1 = dagnode.index;
        int ai[] = new int[2];
        int ai1[] = new int[2];
        float f = (float)i / (float)_xRes;
        float f1 = (float)j / (float)_yRes;
        float f2 = f * 0.5F;
        float f3 = f1 * 0.5F;
        for(int j1 = 0; j1 < dagnode.children.size(); j1++)
        {
            DAGNode dagnode1 = (DAGNode)dagnode.children.elementAt(j1);
            int k1 = dagnode1.index;
            ai[1] = i1 / _xRes;
            ai[0] = i1 - ai[1] * _xRes;
            ai1[1] = k1 / _xRes;
            ai1[0] = k1 - ai1[1] * _xRes;
            ai[1] = _yRes - ai[1];
            ai1[1] = _yRes - ai1[1];
            if(dagnode1.leader)
                graphics2d.setPaint(new Color(leaderCharge, leaderCharge, leaderCharge));
            else
            if(dagnode1.secondary)
            {
                if(dagnode1.maxDepthNode == null)
                {
                    DAGNode dagnode2 = findDeepest(dagnode1);
                    dagnode1.maxDepthNode = dagnode2;
                }
                int l1 = dagnode1.maxDepthNode.depth;
                float f4 = -(float)(l1 * l1) / (float)(Math.log(secondaryCharge) * 2D);
                float f5 = -(float)dagnode1.depth * (float)dagnode1.depth;
                f5 /= 2.0F * f4;
                f5 = (float)Math.exp(f5) * 0.5F;
                graphics2d.setPaint(new Color(f5, f5, f5));
            } else
            {
                graphics2d.setPaint(new Color(0.0F, 1.0F, 0.0F));
            }
            graphics2d.draw(new java.awt.geom.Line2D.Float((float)ai[0] * f + f2 + (float)k, ((float)ai[1] * f1 - f3) + (float)l, (float)ai1[0] * f + f2 + (float)k, ((float)ai1[1] * f1 - f3) + (float)l));
            if(dagnode1.children.size() > 0)
                drawSegments(dagnode1, graphics2d, i, j, k, l);
        }

    }

    public BufferedImage drawTree()
    {
        BackBuffer _tmp = backBuffer;
        char c = '\u0200';
        BackBuffer _tmp1 = backBuffer;
        char c1 = '\u0200';
        BufferedImage bufferedimage = new BufferedImage(c, c1, 1);
        Graphics2D graphics2d;
        try
        {
            graphics2d = (Graphics2D)bufferedimage.getGraphics();
        }
        catch(NullPointerException nullpointerexception)
        {
            nullpointerexception.printStackTrace();
            return null;
        }
        graphics2d.setBackground(Color.black);
        graphics2d.clearRect(0, 0, c, c1);
        drawSegments(root, graphics2d, c, c1, 0, 0);
        backBuffer.backBuffer = bufferedimage;
        backBuffer.backBufferContext = graphics2d;
        return bufferedimage;
    }

    public void run()
    {
        do
        {
            if(pauseFlag)
                pause();
            if(hitBottom())
            {
                pauseFlag = true;
                statusBox.append("Simulation complete.\n");
                statusBox.append("Press 'Add Glow' to apply the final render.\n");
            }
            drawGridAndTree();
            conjugateGradient();
            addParticle();
        } while(true);
    }

    synchronized void wake()
    {
        notify();
    }

    synchronized void pause()
    {
        try
        {
            wait();
        }
        catch(InterruptedException interruptedexception) { 
            interruptedexception.printStackTrace();
        }
    }

    private void conjugateGradient()
    {
        int i5 = 0;
        calcResidual();
        float f = -0.25F;
        int k4 = _xRes + 1;
        for(int j2 = 1; j2 < _yRes - 1;)
        {
            for(int i = 1; i < _xRes - 1;)
            {
                _direction[k4] = _residual[k4] * f;
                i++;
                k4++;
            }

            j2++;
            k4 += 2;
        }

        float f1 = 0.0F;
        k4 = _xRes + 1;
        for(int k2 = 1; k2 < _yRes - 1;)
        {
            for(int j = 1; j < _xRes - 1;)
            {
                f1 += _residual[k4] * _direction[k4];
                j++;
                k4++;
            }

            k2++;
            k4 += 2;
        }

        float f2 = f1;
        float f3 = 1000F;
        float f4 = 0.001F;
        for(float f5 = 1.0F; (float)i5 < f3 && f5 > 1E-05F; f5 = maxResidual())
        {
            int l4 = _xRes + 1;
            for(int l2 = 1; l2 < _yRes - 1;)
            {
                for(int k = 1; k < _xRes - 1;)
                {
                    _q[l4] = -4F * _direction[l4] + _direction[l4 - 1] + _direction[l4 + 1] + _direction[l4 - _xRes] + _direction[l4 + _xRes];
                    k++;
                    l4++;
                }

                l2++;
                l4 += 2;
            }

            float f6 = 0.0F;
            l4 = _xRes + 1;
            for(int i3 = 1; i3 < _yRes - 1;)
            {
                for(int l = 1; l < _xRes - 1;)
                {
                    if(_phase[l4] < 0.5F)
                        f6 += _direction[l4] * _q[l4];
                    l++;
                    l4++;
                }

                i3++;
                l4 += 2;
            }

            f6 = f1 / f6;
            l4 = _xRes + 1;
            for(int j3 = 1; j3 < _yRes - 1;)
            {
                for(int i1 = 1; i1 < _xRes - 1;)
                {
                    if(_phase[l4] < 0.5F)
                        _potential[l4] += f6 * _direction[l4];
                    i1++;
                    l4++;
                }

                j3++;
                l4 += 2;
            }

            l4 = _xRes + 1;
            for(int k3 = 1; k3 < _yRes - 1;)
            {
                for(int j1 = 1; j1 < _xRes - 1;)
                {
                    if(_phase[l4] < 0.5F)
                        _residual[l4] -= f6 * _q[l4];
                    else
                        _residual[l4] = 0.0F;
                    j1++;
                    l4++;
                }

                k3++;
                l4 += 2;
            }

            l4 = _xRes + 1;
            for(int l3 = 1; l3 < _yRes - 1;)
            {
                for(int k1 = 1; k1 < _xRes - 1;)
                {
                    _s[l4] = _residual[l4] * f;
                    k1++;
                    l4++;
                }

                l3++;
                l4 += 2;
            }

            float f7 = f1;
            f1 = 0.0F;
            l4 = _xRes + 1;
            for(int i4 = 1; i4 < _yRes - 1;)
            {
                for(int l1 = 1; l1 < _xRes - 1;)
                {
                    f1 += _residual[l4] * _s[l4];
                    l1++;
                    l4++;
                }

                i4++;
                l4 += 2;
            }

            float f8 = f1 / f7;
            l4 = _xRes + 1;
            for(int j4 = 1; j4 < _yRes - 1;)
            {
                for(int i2 = 1; i2 < _xRes - 1;)
                {
                    if(_phase[l4] < 0.5F)
                        _direction[l4] = _s[l4] + f8 * _direction[l4];
                    else
                        _direction[l4] = 0.0F;
                    i2++;
                    l4++;
                }

                j4++;
                l4 += 2;
            }

            i5++;
        }

    }

    private void calcResidual()
    {
        int k = _xRes + 1;
        for(int j = 1; j < _yRes - 1;)
        {
            for(int i = 1; i < _xRes - 1;)
            {
                if(_phase[k] < 0.5F)
                {
                    _residual[k] = -4F * _potential[k] + _potential[k + 1] + _potential[k - 1] + _potential[k + _xRes] + _potential[k - _xRes];
                    _residual[k] = -_residual[k];
                    _residual[k] += _charge[k] * _dx * _dx;
                } else
                {
                    _residual[k] = 0.0F;
                }
                i++;
                k++;
            }

            j++;
            k += 2;
        }

    }

    private float maxResidual()
    {
        int i = _xRes + 1;
        float f = 0.0F;
        for(int j = 1; j < _yRes - 1;)
        {
            for(int k = 1; k < _xRes - 1;)
            {
                if(_phase[i] < 0.5F)
                {
                    float f1 = -4F * _potential[i] + _potential[i + 1] + _potential[i - 1] + _potential[i + _xRes] + _potential[i - _xRes];
                    f1 = -f1;
                    f1 += _charge[i] * _dx * _dx;
                    f = f1 <= f ? f : f1;
                }
                k++;
                i++;
            }

            j++;
            i += 2;
        }

        return f;
    }

    private void addParticle()
    {
        float f = 1.0F / (float)Math.sqrt(2D);
        Vector vector = new Vector();
        Vector vector1 = new Vector();
        int k = _xRes + 1;
        float f1 = 0.0F;
        for(int j = 1; j < _yRes - 1;)
        {
            for(int i = 1; i < _xRes - 1;)
            {
                boolean flag = _lightning[k + 1] || _lightning[k - 1] || _lightning[k + _xRes] || _lightning[k - _xRes];
                boolean flag1 = _lightning[k + _xRes + 1] || _lightning[(k + _xRes) - 1] || _lightning[k - _xRes - 1] || _lightning[(k - _xRes) + 1];
                if(_phase[k] < 0.5F && flag)
                {
                    vector.add(new Integer(k));
                    float f3 = (float)Math.pow(_potential[k], _power);
                    f1 += f3;
                    vector1.add(new Float(f3));
                } else
                if(_phase[k] < 0.5F && flag1)
                {
                    vector.add(new Integer(k));
                    float f4 = (float)Math.pow(_potential[k] * f, _power);
                    f1 += f4;
                    vector1.add(new Float(f4));
                }
                i++;
                k++;
            }

            j++;
            k += 2;
        }

        if(vector.size() == 0)
            return;
        int l = 0;
        if(f1 < (float)Math.pow(1.000000013351432E-10D, _power))
        {
            l = (int)((float)vector.size() * rand.nextFloat());
        } else
        {
            float f2 = rand.nextFloat();
            float f5 = 1.0F / f1;
            float f6 = ((Float)vector1.elementAt(0)).floatValue();
            for(float f7 = f6 * f5; f7 < f2 && l < vector.size(); f7 += f6 * f5)
            {
                l++;
                f6 = ((Float)vector1.elementAt(l)).floatValue();
            }

        }
        int i1 = ((Integer)vector.elementAt(l)).intValue();
        _phase[i1] = 1.0F;
        _potential[i1] = 0.0F;
        _lightning[i1] = true;
        addSegment(i1);
    }

    private boolean hitBottom()
    {
        for(int i = _xRes + 1; i < _xRes + (_xRes - 1); i++)
        {
            int j = i;
            if(_lightning[j])
            {
                bottomHit = j;
                buildLeader();
                return true;
            }
        }

        return false;
    }

    private void buildLeader()
    {
        DAGNode dagnode = (DAGNode)DAGHash.get(new Integer(bottomHit));
        for(int i = 0; i < _xRes * _yRes; i++)
            _charge[i] = 0.0F;

        for(; dagnode != null; dagnode = dagnode.parent)
        {
            dagnode.leader = true;
            dagnode.secondary = false;
            for(int j = 0; j < dagnode.children.size(); j++)
            {
                DAGNode dagnode1 = (DAGNode)dagnode.children.elementAt(j);
                if(dagnode1.secondary)
                {
                    buildBranch(dagnode1, 1);
                    buildSecondaryLongest(dagnode1);
                }
            }

        }

    }

    private void buildBranch(DAGNode dagnode, int i)
    {
        dagnode.depth = i;
        dagnode.leader = false;
        for(int j = 0; j < dagnode.children.size(); j++)
        {
            DAGNode dagnode1 = (DAGNode)dagnode.children.elementAt(j);
            if(dagnode1.secondary)
                buildBranch(dagnode1, i + 1);
        }

    }

    public void buildSecondaryLongest(DAGNode dagnode)
    {
        DAGNode dagnode1 = findDeepest(dagnode);
        if((float)dagnode1.depth * _dx < longestBranch)
            return;
        for(DAGNode dagnode2 = dagnode1; dagnode2 != dagnode; dagnode2 = dagnode2.parent)
        {
            dagnode2.leader = false;
            dagnode2.secondary = true;
            buildNtharyLongest(dagnode2);
        }

        dagnode.secondary = true;
    }

    public void buildNtharyLongest(DAGNode dagnode)
    {
        for(int i = 0; i < dagnode.children.size(); i++)
        {
            DAGNode dagnode1 = (DAGNode)dagnode.children.elementAt(i);
            if(dagnode1.secondary)
                continue;
            DAGNode dagnode2 = findDeepest(dagnode1);
            int j = dagnode2.depth - dagnode.depth;
            if(j <= 20)
                continue;
            for(DAGNode dagnode3 = dagnode2; dagnode3 != dagnode; dagnode3 = dagnode3.parent)
                dagnode3.secondary = true;

            dagnode.secondary = true;
        }

    }

    public DAGNode findDeepest(DAGNode dagnode)
    {
        DAGNode dagnode1 = dagnode;
        for(int i = 0; i < dagnode.children.size(); i++)
        {
            DAGNode dagnode2 = (DAGNode)dagnode.children.elementAt(i);
            DAGNode dagnode3 = findDeepest(dagnode2);
            if(dagnode3.depth > dagnode1.depth)
                dagnode1 = dagnode3;
        }

        return dagnode1;
    }

    public void setCharge(int i, int j)
    {
        int k = i + j * _xRes;
        _charge[k] = 2.0F;
        _phase[k] = 1.0F;
        _potential[k] = 1.0F;
        drawGridAndTree();
    }

    public boolean pauseFlag;
    private int _xRes;
    private int _yRes;
    private int _totalCells;
    private float _width;
    private float _height;
    private float _dx;
    private float _dy;
    private float _phase[];
    private float _potential[];
    private float _charge[];
    private boolean _lightning[];
    private float _power;
    private float _residual[];
    private float _direction[];
    private float _q[];
    private float _s[];
    BackBuffer backBuffer;
    TextArea statusBox;
    Random rand;
    DAGNode root;
    HashMap DAGHash;
    Vector nodeOrder;
    int bottomHit;
    APSF apsf;
    int apsfRes;
    float leaderCharge;
    float secondaryCharge;
    float sideCharge;
    float longestBranch;
    public int chargeType;
}

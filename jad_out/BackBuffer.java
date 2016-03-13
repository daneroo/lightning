// Decompiled by Jad v1.5.8e. Copyright 2001 Pavel Kouznetsov.
// Jad home page: http://www.geocities.com/kpdus/jad.html
// Decompiler options: packimports(3) 
// Source File Name:   BackBuffer.java

import java.awt.*;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.image.BufferedImage;

public class BackBuffer extends Canvas
    implements Runnable, MouseListener
{

    public void setLightning(Lightning lightning1)
    {
        lightning = lightning1;
    }

    BackBuffer()
    {
        bgColor = new Color(0, 0, 100);
        backBuffer = null;
        backBufferContext = null;
        lightning = null;
        pauseFlag = false;
        setSize(512, 512);
    }

    public void run()
    {
        if(backBufferContext == null)
            init();
        Graphics g = getGraphics();
        do
        {
            paintBuffer(g);
            if(pauseFlag)
                pause();
            try
            {
                Thread.sleep(10L);
            }
            catch(InterruptedException interruptedexception) { 
                interruptedexception.printStackTrace();
            }
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

    public void paint(Graphics g)
    {
        if(backBuffer != null)
            g.drawImage(backBuffer, 0, 0, this);
    }

    public void update(Graphics g)
    {
        paint(g);
    }

    public void init()
    {
        setSize(512, 512);
        addMouseListener(this);
        try
        {
            backBuffer = new BufferedImage(512, 512, 1);
            backBufferContext = backBuffer.getGraphics();
        }
        catch(Exception exception)
        {
            backBuffer = null;
            backBufferContext = null;
            return;
        }
    }

    void paintBuffer(Graphics g)
    {
        g.drawImage(backBuffer, 0, 0, this);
    }

    public void mouseClicked(MouseEvent mouseevent)
    {
        if(lightning.chargeType == -1)
            return;
        int i = mouseevent.getX();
        int j = mouseevent.getY();
        if(i < 512 && i > 0 && j < 512 && j > 0)
        {
            int k = (int)(((float)i / 512F) * (float)lightning.xRes());
            int l = (int)(((float)j / 512F) * (float)lightning.yRes()) + 1;
            lightning.setCharge(k, lightning.yRes() - l);
        }
        mouseevent.consume();
    }

    public void mouseReleased(MouseEvent mouseevent)
    {
    }

    public void mousePressed(MouseEvent mouseevent)
    {
    }

    public void mouseExited(MouseEvent mouseevent)
    {
    }

    public void mouseEntered(MouseEvent mouseevent)
    {
    }

    static final int CANVASWIDTH = 512;
    static final int CANVASHEIGHT = 512;
    Color bgColor;
    public boolean pauseFlag;
    BufferedImage backBuffer;
    Graphics backBufferContext;
    Lightning lightning;
}

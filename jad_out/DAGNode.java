// Decompiled by Jad v1.5.8e. Copyright 2001 Pavel Kouznetsov.
// Jad home page: http://www.geocities.com/kpdus/jad.html
// Decompiler options: packimports(3) 
// Source File Name:   Lightning.java

import java.util.Vector;

class DAGNode
{

    DAGNode(int i)
    {
        index = i;
        parent = null;
        leader = true;
        secondary = true;
        depth = 0;
        children = new Vector();
        maxDepthNode = null;
    }

    int index;
    Vector children;
    DAGNode parent;
    boolean leader;
    boolean secondary;
    int depth;
    float potential;
    DAGNode maxDepthNode;
}

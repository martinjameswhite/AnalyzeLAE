#!/usr/bin/env python
import numpy as np


def rotate_to(pra,pdc,cra,cdc):
    """Rotates an array of points -- (RA,DEC) in degrees -- such
    that (0,0) is mapped to (CRA,CDC), also in degrees."""
    tt,pp= np.radians(90-pdc),np.radians(pra)
    nhat = np.zeros( (3,tt.size) )
    nhat[0,:] = np.sin(tt)*np.cos(pp)
    nhat[1,:] = np.sin(tt)*np.sin(pp)
    nhat[2,:] = np.cos(tt)
    # Normally we would rotate the north pole to the final
    # position, but (RA,DEC)=(0,0) is actually the x-axis.
    # This means we set theta->90-theta.
    tt,pp= np.radians(cdc),np.radians(cra)
    rot1 = np.array([[np.cos(tt),0,-np.sin(tt)],\
                     [    0     ,1,     0     ],\
                     [np.sin(tt),0,np.cos(tt)]])
    rot2 = np.array([[np.cos(pp),-np.sin(pp),0],\
                     [np.sin(pp), np.cos(pp),0],\
                     [   0      ,     0     ,1]])
    #
    nhat = np.dot(rot2,np.dot(rot1,nhat))
    rdc  = 90 - np.arccos(nhat[2,:]) * 180./np.pi
    rra  = np.arctan2(nhat[1,:],nhat[0,:]) * 180./np.pi
    return( (rra,rdc) )

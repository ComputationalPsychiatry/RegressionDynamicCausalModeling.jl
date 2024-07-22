/*%% --------------------------------------------------------------------------------------------------
% dcm_euler_integration - Integrates the DCM-fMRI dynamical system using Euler's method
% [x,s,f,v,q]  = dcm_euler_integration(A,C,U,B,D,...
                   rho,alphaInv,tau,gamma,kappa,paramList);

%---------------------------------------------------------------------------------------------
% INPUT:
 * Pls note: All matrices within A, B, D are transposes of the original matrix
%       A           - Region connection matrix
%       C           - Input contribution, represents U*C' (for optimization purposes)
%       U           - Input matrix
%       B           - Bi linear connection matrix
%       D           - Non linear connection matrix
%       rho         - hemodynamic const - one for each region (array of nStates)
%       alphaInv    - hemodynamic const
%       tau         - one for each region (array of nStates)
%       gamma       - hemodynamic const
%       kappa       - one for each region (array of nStates)
%       paramList   - [timeStep nTime nStates nInputs subjectNo] 
 *                      timeStep - euler step size
 *                      nTime - total steps
 *                      nStates - number of regions
 *                      nInputs - number of inputs
 *                      subjectNo- not used by function
%
% Optional:
%
%--------------------------------------------------------------------------------------------
% OUTPUT:
 *      x - neural activity
 *      s - vasodilatory signal
 *      f - blood flow
 *      v - blood volume
 *      q - deoxyhemoglobin content
%           
% -------------------------------------------------------------------------------------------
% REFERENCE:
%
% Author:   Sudhir Shankar Raman, TNU, UZH & ETHZ - 2013
% Modified: Imre Kertesz, TNU, UZH & ETHZ - 2023
%
% This file is part of the RegressionDynamicCausalModeling.jl toolbox, which is part of
% TAPAS, a collection of software tools to support clinical neuromodeling, particularly
% computational psychiatry, computational neurology, and computational psychosomatics.

%%*/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>


void dcm_euler_integration(double *A, double *C, double *U, double *B, double *D, double *rho,
                        double alphaInv, double *tau, double gamma, double *kappa, double *param,
                        double *x_out,double *s_out,double *f1_out,double *v1_out,double *q1_out)
{

    int jState,mState,kIter;
    long iStep;
    double temp1; 
    double temp2;
    double timeStep = param[0];   // time step (U.u.dt)
    long nTime = param[1];        // size(U,1)
    const int nStates = param[2]; // number of regions
    int nInputs = param[3];       // number of inputs
    int dcmTypeB = param[5];      // whether or not to calculate bi-linear part
    int dcmTypeD = param[6];      // if we should also calculate non-linearities
    double *oldf1 = (double *) malloc(nStates*sizeof(double));
    double *oldv1 = (double *) malloc(nStates*sizeof(double));
    double *oldq1 = (double *) malloc(nStates*sizeof(double));
    
    /* Initialize the dynamical system to resting state values */
    for(jState=0;jState<nStates;jState++)
    {
        x_out[nTime*jState] = 0.0;
        s_out[nTime*jState] = 0.0;
        f1_out[nTime*jState] = 1.0;
        v1_out[nTime*jState] = 1.0;
        q1_out[nTime*jState] = 1.0;          
        oldf1[jState] = 0.0;
        oldv1[jState] = 0.0;
        oldq1[jState] = 0.0;
    }
    
    /* Euler's integration steps*/
    for(iStep=0;iStep<(nTime-1);iStep++)
    /*for(iStep=0;iStep<25;iStep++)*/
    {        
        /* For each region*/
        for(jState=0;jState<nStates;jState++)
        {
            /* update x */
            x_out[iStep+1 + nTime*jState] = 0.0;
            temp1 = 0.0;            
            for(kIter=0;kIter<nStates;kIter++)
            {
                temp1=temp1+x_out[iStep + nTime*kIter]*A[kIter + nStates*jState];
            }            
            x_out[iStep+1 + nTime*jState]= x_out[iStep + nTime*jState] + timeStep* (temp1+ C[iStep + nTime*jState]); 
            
            if(dcmTypeB != 0)
            {
                /*B matrix update*/
                for(kIter=0;kIter<nInputs;kIter++)
                {
                    temp1 = 0.0;
                    for(mState=0;mState<nStates;mState++)
                    {
                        temp1=temp1+x_out[iStep + nTime*mState]*B[mState + nStates*(jState + nStates*kIter)];
                    }
                    x_out[iStep+1 + nTime*jState] = x_out[iStep+1 + nTime*jState] + timeStep*(U[iStep +nTime*kIter]*temp1);
                }
            }
            
            if(dcmTypeD != 0)
            {
                /*D matrix update*/
                for(kIter=0;kIter<nStates;kIter++)         
                {
                     temp1 = 0.0;
                     for(mState=0;mState<nStates;mState++)
                     {
                         temp1=temp1+x_out[iStep + nTime*mState]*D[mState + nStates*(jState + nStates*kIter)];
                     }
                     x_out[iStep+1 + nTime*jState] = x_out[iStep+1 + nTime*jState] + timeStep*(x_out[iStep +nTime*kIter]*temp1);
                }
            }
            /* update s */
            s_out[iStep+1 + nTime*jState] = s_out[iStep + nTime*jState] + timeStep*(x_out[iStep + nTime*jState] -
                                            kappa[jState]*s_out[iStep + nTime*jState] - 
                                            gamma*(f1_out[iStep + nTime*jState] - 1.0));
            
            /*if((iStep>15) && (iStep <28))
            {
                printf("%e  %e  %e %e\n",gamma*(f1_out[iStep + nTime*jState] - 1.0),kappa[jState]*s_out[iStep + nTime*jState],x_out[iStep + nTime*jState],s_out[iStep+1 + nTime*jState]);
            }*/
            
            /* update f */
            oldf1[jState]= oldf1[jState] + timeStep*(s_out[iStep + nTime*jState]/f1_out[iStep + nTime*jState]);
            /*f1_out[iStep+1 + nTime*jState]= f1_out[iStep + nTime*jState] + timeStep*(s_out[iStep + nTime*jState]);*/
            
            /* update v */
           temp1 = pow(v1_out[iStep + nTime*jState],alphaInv);
           temp2 = (1.0- pow((1.0-rho[jState]),(1.0/f1_out[iStep + nTime*jState])))/rho[jState];
            oldv1[jState]= oldv1[jState] + timeStep*((f1_out[iStep + nTime*jState] - temp1)/
                           (tau[jState]*v1_out[iStep + nTime*jState]));
            
            /* update q */
            oldq1[jState]= oldq1[jState] + timeStep*((f1_out[iStep + nTime*jState]*temp2 - 
                    temp1*q1_out[iStep + nTime*jState]/v1_out[iStep + nTime*jState])/(tau[jState]*q1_out[iStep + nTime*jState]));

            /* tracking the exponentiated values */
            f1_out[iStep+1 + nTime*jState] = exp(oldf1[jState]);
            v1_out[iStep+1 + nTime*jState] = exp(oldv1[jState]);
            q1_out[iStep+1 + nTime*jState] = exp(oldq1[jState]);
        }       
    }
    
    /* free all the created memory*/
    free(oldf1);
    free(oldv1);
    free(oldq1);
    return;    
}
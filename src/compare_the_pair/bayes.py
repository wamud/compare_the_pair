# Some stuff Alan dumped here
import random
import math
import numpy as np
import copy


import sys, os, time

from qinfer import LiuWestResampler
from qinfer import utils


class Distribution():
    
    def __init__(self, distribution_generator=None,
                 n_points=100,
                 n_qubits=1, 
                 rejection_threshold=0.5,
                 resampler_a=np.sqrt(0.9),
                 ):
        
        
        if distribution_generator is not None:
            self.points, self.weights = distribution_generator(n_points, n_qubits)
        else:
            self.points, self.weights = self.random_distribution(n_points=n_points)
        
        self.resampler_a = resampler_a
        
        self.resampler = LiuWestResampler(a=self.resampler_a, postselect=False)
        self.n_qubits = n_qubits
        self.n_points = len(self.weights)
        self.rejection_threshold = rejection_threshold
        self.covariance_matrix = None
        
        
        
    def update(self, measurement_data, resample=True):
        self.update_estimate(measurement_data)
        return

    def renormalise(self):
        self.weights /= sum(self.weights)
        
        
    def calc_covariance_matrix(self):
        ''' Calculates the covariance matrix
        '''
        return utils.particle_covariance_mtx(self.weights, self.points)     
        
    
    def resample(self):

        self.renormalise()
        self.points = self.points.reshape(len(self.points), 1)
        self.weights, self.points = self.resampler(None, self.weights, self.points)

        self.points = self.points.reshape(self.points.shape[0]) 
        
        
    def update_estimate(self, measurement_outcomes):
        '''
            Calculate the new weight from the previous one
        '''        
        for state in enumerate(measurement_outcomes):
            self.conditional_resample()
            if state == 1:
                self.weights *= self.points 
            else:
                self.weights *= (1 - self.points) 
        self.renormalise() 

    def calc_bayes_mean(self):
        return np.sum(self.points * self.weights)
            
        
    def random_distribution(self, n_points, start=0, stop=0.05):

        points = (np.random.rand(n_points) * (stop - start) + start)
        
            
        weights = np.array([1 / len(points)] * len(points))
        points = np.array(points)
        
        return points, weights
        
    
    def calc_bayes_risk(self):
        '''
        Simple model using mean squared error as the loss function
        This sets the estimate to the bayes mean of the posterior
        '''
        estimate = self.calc_bayes_mean()
        return np.sum(self.weights * np.abs(self.points - estimate) ** 2)
        
    
    def n_eff(self):
        return 1 / sum(i ** 2 for i in self.weights)
    
    def conditional_resample(self):
        
        # Calculate number of effective particles        
        if self.n_eff() < self.n_points * self.rejection_threshold:
            self.resample()
    
    
        
    def linear_distribution(self, step, start=0, stop=1, n_qubits=1):
        '''
            Simple linear, evenly weighted distribution
        '''
        points = []
        
        position = start
        while position <= stop:
            points.append(position)
            position += step
        
        if stop > points[-1]:
            points.append(stop)
            
        weights = np.array([1 / len(points)] * len(points))
        
        points = np.array(points)
        
        return points, weights


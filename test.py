import tensorflow as tf
import numpy as np
import math

############################################################################################################

# Import csv file, then generate numpy array
from numpy import genfromtxt
data = genfromtxt('/Users/JYang/Desktop/Stanford/gene_effects_data.csv', delimiter=',')

# Python optimisation variables
learning_rate = 0.5 #probably do 0.01 or 0.1 or 0.05
epochs = 10 #check size of final dataset
batch_size = 100 #check size of final dataset

#get shape of file
dim = data.shape

# declare the training data placeholders
# input x - number of genes (? x 52408) --> (#samples x #features)
x = tf.placeholder(tf.float32, [None, dim[0]])
# now declare the output data placeholder - binary output (0, 1) --> (protective, susceptible)
y = tf.placeholder(tf.float32, [None, 2])

#weights connecting input layer to hidden layer
#intilized with random normal distribution (mean of 0, std dev of ...?)
#start by testing a hidden layer that's ~half the size of the input
W1 = tf.Variable(tf.random_normal([dim[1], int(math.ceil(dim[1]/2))], stddev = 0.03), name = 'W1')
b1 = tf.Variable(tf.random_normal([int(math.ceil(dim[1]/2))], name = 'b1'))
#TODO: add another hidden layer, depending on ability of 3 layer neural net
#weights connecting hidden layer to output layer
W2 = tf.Variable(tf.random_normal([int(math.ceil(dim[1]/2)), 2], stddev = 0.03), name = 'W2')
b2 = tf.Variable(tf.random_normal([10]), name='b2')

#setup input nodes
#define activation functions for hidden layers (rectified linear unit, sigmoid, tanh...)
hidden_out = tf.add(tf.matmul(x, W1), b1)
#try rectified hidden linear unit
hidden_out = tf.nn.relu(hidden_out)

#get the output layer
y_ = tf.add(tf.matmul(hidden_out, W2), b2)
y_ = tf.nn.softmax(y_)

#clip output so we never have log(0)
y_clipped = tf.clip_by_value(y_, 1e-10, 0.9999999)
#loss function (cross entropy, mean squared error...)
cross_entropy = -tf.reduce_mean(tf.reduce_sum(y * tf.log(y_clipped) + (1 - y) * tf.log(1 - y_clipped), axis=1))

#set up an optimizer
optimizer = tf.train.GradientDescentOptimizer(learning_rate = learning_rate).minimize(cross_entropy)

#set up variable initialization
init_op = tf.global_variables_initializer()

#define accuracy assessment operation
correct_prediction = tf.equal(tf.argmax(y, 1), tf.argmax(y_, 1))
accuracy = tf.reduce_mean(tf.cast(correct_prediction, tf.float32))

#set up and start the session
with tf.Session() as sess:
   # initialise the variables
   sess.run(init_op)
   #edit from here on
   total_batch = int(len(mnist.train.labels) / batch_size)
   for epoch in range(epochs):
        avg_cost = 0
        for i in range(total_batch):
            batch_x, batch_y = mnist.train.next_batch(batch_size=batch_size)
             _, c = sess.run([optimiser, cross_entropy], 
                         feed_dict={x: batch_x, y: batch_y})
            avg_cost += c / total_batch
        print("Epoch:", (epoch + 1), "cost =", "{:.3f}".format(avg_cost))
   print(sess.run(accuracy, feed_dict={x: mnist.test.images, y: mnist.test.labels}))

#TODO: def load_data()

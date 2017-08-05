import tensorflow as tf

from tensorflow.examples.tutorials.mnist import input_data
mnist = input_data.read_data_sets('MNIST_data', one_hot=True)
import numpy as np
x = tf.placeholder(tf.float32, [None, 784])
W = tf.Variable(tf.zeros([784, 10]))
b = tf.Variable(tf.zeros([10]))
y = tf.nn.softmax(tf.matmul(x, W) + b)
y_ = tf.placeholder(tf.float32, [None, 10])
cross_entropy = tf.reduce_mean(-tf.reduce_sum(y_ * tf.log(y), reduction_indices=[1]))
train_step = tf.train.GradientDescentOptimizer(0.5).minimize(cross_entropy)
init = tf.global_variables_initializer()
sess = tf.Session()
sess.run(init)
for i in range(1000):
  batch_xs, batch_ys = mnist.train.next_batch(100)
  sess.run(train_step, feed_dict={x: batch_xs, y_: batch_ys})
  if(i%100==0):
    correct_prediction = tf.equal(tf.argmax(y,1), tf.argmax(y_,1))
    accuracy = tf.reduce_mean(tf.cast(correct_prediction, tf.float32))
    print(sess.run(accuracy, feed_dict={x: mnist.test.images, y_: mnist.test.labels}))

'''
fileName = '../../../data/seis/dfault/fake/2d/'
txs = np.fromfile(fileName+"data61.dat", dtype=np.single)
tys = np.fromfile(fileName+"mark61.dat", dtype=np.single)
txs = txs.reshape(100,600)
tys = tys.reshape(100,22)
for i in range(600):
  xs = np.fromfile(fileName+"data"+str(i)+".dat", dtype=np.single)
  ys = np.fromfile(fileName+"mark"+str(i)+".dat", dtype=np.single)
  xs = xs.reshape(100,600)
  ys = ys.reshape(100,22)
  sess.run(train_step, feed_dict={x: xs, y_: ys})
  if(i%50==0):
    correct_prediction = tf.equal(tf.argmax(y,1), tf.argmax(y_,1))
    accuracy = tf.reduce_mean(tf.cast(correct_prediction, tf.float32))
    print(sess.run(accuracy, feed_dict={x: txs, y_: tys}))
'''

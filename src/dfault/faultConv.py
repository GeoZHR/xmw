import tensorflow as tf
import math 
import struct

#from tensorflow.examples.tutorials.mnist import input_data
#mnist = input_data.read_data_sets('MNIST_data', one_hot=True, reshape=False)
import numpy as np

fileName = '../../../data/seis/dfault/fake/2d/'
txs = np.fromfile(fileName+"seis1.dat", dtype=np.single)
tys = np.fromfile(fileName+"marks1.dat", dtype=np.single)
txs = txs.reshape(100,28,28)
tys = tys.reshape(100,22)
txsa = np.ndarray(shape=(100,28,28,1),dtype=np.single)
n1,n2,n3=28,28,100
for i3 in range(n3):
  for i2 in range(n2):
    for i1 in range(n1):
      txsa[i3][i2][i1][0] = txs[i3][i2][i1]
# input X: 28x28 grayscale images, the first dimension (None) will index the images in the mini-batch
X = tf.placeholder(tf.float32, [None, 28, 28, 1])
# correct answers will go here
Y_ = tf.placeholder(tf.float32, [None, 22])
# variable learning rate
lr = tf.placeholder(tf.float32)
# test flag for batch norm
iter = tf.placeholder(tf.int32)
tst = tf.placeholder(tf.bool)

# dropout probability
pkeep_conv = tf.placeholder(tf.float32)
pkeep = tf.placeholder(tf.float32)

def batchnorm(Ylogits, is_test, iteration, offset, convolutional=False):
    exp_moving_avg = tf.train.ExponentialMovingAverage(0.999, iteration) # adding the iteration prevents from averaging across non-existing iterations
    bnepsilon = 1e-5
    if convolutional:
        mean, variance = tf.nn.moments(Ylogits, [0, 1, 2])
    else:
        mean, variance = tf.nn.moments(Ylogits, [0])
    update_moving_everages = exp_moving_avg.apply([mean, variance])
    m = tf.cond(is_test, lambda: exp_moving_avg.average(mean), lambda: mean)
    v = tf.cond(is_test, lambda: exp_moving_avg.average(variance), lambda: variance)
    Ybn = tf.nn.batch_normalization(Ylogits, m, v, offset, None, bnepsilon)
    return Ybn, update_moving_everages

def no_batchnorm(Ylogits, is_test, iteration, offset, convolutional=False):
    return Ylogits, tf.no_op()

def compatible_convolutional_noise_shape(Y):
    noiseshape = tf.shape(Y)
    noiseshape = noiseshape * tf.constant([1,0,0,1]) + tf.constant([0,1,1,0])
    return noiseshape


# three convolutional layers with their channel counts, and a
# fully connected layer (tha last layer has 10 softmax neurons)
L1 = 64  # first convolutional layer output depth
L2 = 64  # second convolutional layer output depth
L3 = 64  # third convolutional layer
L4 = 128  # fully connected layer
L5 = 128  # fully connected layer
L6 = 138  # fully connected layer
L7 = 256  # fully connected layer
L8 = 256  # fully connected layer


W1 = tf.Variable(tf.truncated_normal([3, 3, 1, L1], stddev=0.1))  # 6x6 patch, 1 input channel, K output channels
B1 = tf.Variable(tf.constant(0.1, tf.float32, [L1]))

W2 = tf.Variable(tf.truncated_normal([3, 3, L1, L2], stddev=0.1))
B2 = tf.Variable(tf.constant(0.1, tf.float32, [L2]))

W3 = tf.Variable(tf.truncated_normal([3, 3, L2, L3], stddev=0.1))
B3 = tf.Variable(tf.constant(0.1, tf.float32, [L3]))

W4 = tf.Variable(tf.truncated_normal([3, 3, L3, L4], stddev=0.1))
B4 = tf.Variable(tf.constant(0.1, tf.float32, [L4]))

W5 = tf.Variable(tf.truncated_normal([3, 3, L4, L5], stddev=0.1))
B5 = tf.Variable(tf.constant(0.1, tf.float32, [L5]))

W6 = tf.Variable(tf.truncated_normal([3, 3, L5, L6], stddev=0.1))
B6 = tf.Variable(tf.constant(0.1, tf.float32, [L6]))

W7 = tf.Variable(tf.truncated_normal([3, 3, L6, L7], stddev=0.1))
B7 = tf.Variable(tf.constant(0.1, tf.float32, [L7]))

W8 = tf.Variable(tf.truncated_normal([7 * 7 * L7, L8], stddev=0.1))
B8 = tf.Variable(tf.constant(0.1, tf.float32, [L8]))

W9 = tf.Variable(tf.truncated_normal([L8, 22], stddev=0.1))
B9 = tf.Variable(tf.constant(0.1, tf.float32, [22]))


# The model
# batch norm scaling is not useful with relus
# batch norm offsets are used instead of biases
pooling_kernel=(1,2,2,1)
stride = 1  # output is 28x28
Y1l = tf.nn.conv2d(X, W1, strides=[1, stride, stride, 1], padding='SAME')
Y1bn, update_ema1 = batchnorm(Y1l, tst, iter, B1, convolutional=True)
Y1 = tf.nn.relu(Y1bn)
#Y1 = tf.nn.dropout(Y1r, pkeep_conv, compatible_convolutional_noise_shape(Y1r))

stride = 1  # output is 14x14
Y2l = tf.nn.conv2d(Y1, W2, strides=[1, stride, stride, 1], padding='SAME')
Y2bn, update_ema2 = batchnorm(Y2l, tst, iter, B2, convolutional=True)
Y2 = tf.nn.relu(Y2bn)
#Y2 = tf.nn.dropout(Y2r, pkeep_conv, compatible_convolutional_noise_shape(Y2r))

stride = 1  # output is 7x7
Y3l = tf.nn.conv2d(Y2, W3, strides=[1, stride, stride, 1], padding='SAME')
Y3bn, update_ema3 = batchnorm(Y3l, tst, iter, B3, convolutional=True)
Y3r = tf.nn.relu(Y3bn)
#Y3 = tf.nn.dropout(Y3r, pkeep_conv, compatible_convolutional_noise_shape(Y3r))
Y3 = tf.nn.max_pool(Y3r, ksize=pooling_kernel, strides=pooling_kernel,padding='SAME')


stride = 1  # output is 7x7
Y4l = tf.nn.conv2d(Y3, W4, strides=[1, stride, stride, 1], padding='SAME')
Y4bn, update_ema4 = batchnorm(Y4l, tst, iter, B4, convolutional=True)
Y4 = tf.nn.relu(Y4bn)
#Y4 = tf.nn.dropout(Y4r, pkeep_conv, compatible_convolutional_noise_shape(Y4r))

stride = 1  # output is 7x7
Y5l = tf.nn.conv2d(Y4, W5, strides=[1, stride, stride, 1], padding='SAME')
Y5bn, update_ema5 = batchnorm(Y5l, tst, iter, B5, convolutional=True)
Y5 = tf.nn.relu(Y5bn)
#Y5 = tf.nn.dropout(Y5r, pkeep_conv, compatible_convolutional_noise_shape(Y5r))

stride = 1  # output is 7x7
Y6l = tf.nn.conv2d(Y5, W6, strides=[1, stride, stride, 1], padding='SAME')
Y6bn, update_ema6 = batchnorm(Y6l, tst, iter, B6, convolutional=True)
Y6 = tf.nn.relu(Y6bn)
#Y6 = tf.nn.dropout(Y6r, pkeep_conv, compatible_convolutional_noise_shape(Y6r))

stride = 1  # output is 7x7
Y7l = tf.nn.conv2d(Y6, W7, strides=[1, stride, stride, 1], padding='SAME')
Y7bn, update_ema7 = batchnorm(Y7l, tst, iter, B7, convolutional=True)
Y7r = tf.nn.relu(Y7bn)
#Y7 = tf.nn.dropout(Y7r, pkeep_conv, compatible_convolutional_noise_shape(Y7r))
Y7 = tf.nn.max_pool(Y7r, ksize=pooling_kernel, strides=pooling_kernel,padding='SAME')

# reshape the output from the third convolution for the fully connected layer
YY = tf.reshape(Y7, shape=[-1, 7 * 7 * L7])

Y8l = tf.matmul(YY, W8)
Y8bn, update_ema8 = batchnorm(Y8l, tst, iter, B8)
Y8 = tf.nn.relu(Y8bn)
#Y8 = tf.nn.dropout(Y8r, pkeep)
Ylogits = tf.matmul(Y8, W9) + B9
Y = tf.nn.softmax(Ylogits)

update_ema = tf.group(update_ema1, update_ema2, update_ema3, update_ema4)

# cross-entropy loss function (= -sum(Y_i * log(Yi)) ), normalised for batches of 100  images
# TensorFlow provides the softmax_cross_entropy_with_logits function to avoid numerical stability
# problems with log(0) which is NaN
cross_entropy = tf.nn.softmax_cross_entropy_with_logits(logits=Ylogits, labels=Y_)
cross_entropy = tf.reduce_mean(cross_entropy)*100

# accuracy of the trained model, between 0 (worst) and 1 (best)
predict = Y
correct_prediction = tf.equal(tf.argmax(Y, 1), tf.argmax(Y_, 1))
accuracy = tf.reduce_mean(tf.cast(correct_prediction, tf.float32))


# training step, the learning rate is a placeholder
train_step = tf.train.AdamOptimizer(lr).minimize(cross_entropy)

# init
init = tf.global_variables_initializer()
sess = tf.Session()
sess.run(init)

for i in range(1000):
  max_learning_rate = 0.003
  min_learning_rate = 0.0001
  decay_speed = 1600.0
  learning_rate = min_learning_rate + (max_learning_rate - min_learning_rate) * math.exp(-i/decay_speed)
  xs = np.fromfile(fileName+"seis"+str(i)+".dat", dtype=np.single)
  ys = np.fromfile(fileName+"marks"+str(i)+".dat", dtype=np.single)
  ys = ys.reshape(100,22)
  xs = xs.reshape(100,28,28)
  xsa = np.ndarray(shape=(100,28,28,1),dtype=np.single)
  n1,n2,n3=28,28,100
  for i3 in range(n3):
    for i2 in range(n2):
      for i1 in range(n1):
        xsa[i3][i2][i1][0] = xs[i3][i2][i1]
  #batch_xs, batch_ys = mnist.train.next_batch(100)
  #sess.run(train_step, {X: batch_xs, Y_: batch_ys, lr: learning_rate})
  sess.run(train_step, {X: xsa, Y_: ys, lr: learning_rate, tst: False, pkeep: 0.75, pkeep_conv: 1.0})
  sess.run(update_ema, {X: xsa, Y_: ys, tst: False, iter: i, pkeep: 1.0, pkeep_conv: 1.0})
  #sess.run(train_step, {X: xsa, Y_: ys, lr: learning_rate})
  if(i%20==0):
    #a, c = sess.run([accuracy, cross_entropy], {X: batch_xs, Y_: batch_ys})
    #a, c = sess.run([accuracy, cross_entropy], {X: xsa, Y_: ys})ko
    p,a,c = sess.run([predict, accuracy, cross_entropy], {X: txsa, Y_: tys, tst: True, iter: i, pkeep: 1.0, pkeep_conv: 1.0})
    p.tofile(fileName+"label"+str(i)+".dat", format="%4")
    print(str(i) + ": accuracy:" + str(a) + " loss: " + str(c))

'''
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

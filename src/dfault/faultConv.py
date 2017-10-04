# encoding: UTF-8
# Copyright 2016 Google.com
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import tensorflow as tf
import math
import numpy as np

fileName = '../../../data/seis/dfault/fake/2d/'
txs = np.fromfile(fileName+"seisd1.dat", dtype=np.single)
tys = np.fromfile(fileName+"labeld1.dat", dtype=np.single)
txs = txs.reshape(100,32,32)
tys = tys.reshape(100,23)
txsa = np.ndarray(shape=(100,32,32,1),dtype=np.single)
n1,n2,n3=32,32,100
for i3 in range(n3):
  for i2 in range(n2):
    for i1 in range(n1):
      txsa[i3][i2][i1][0] = txs[i3][i2][i1]

#seist = np.fromfile(fileName+"seist50.dat", dtype=np.single)
#seist = seist.reshape(10302,32,32)
#seista = np.ndarray(shape=(10302,32,32,1),dtype=np.single)
#n1,n2,n3=32,32,10302
seist = np.fromfile(fileName+"seisF3d.dat", dtype=np.single)
seist = seist.reshape(818,32,32)
seista = np.ndarray(shape=(818,32,32,1),dtype=np.single)
n1,n2,n3=32,32,818
for i3 in range(n3):
  for i2 in range(n2):
    for i1 in range(n1):
      seista[i3][i2][i1][0] = seist[i3][i2][i1]
# neural network structure for this sample:
#
# · · · · · · · · · ·      (input data, 1-deep)                 X [batch, 28, 28, 1]
# @ @ @ @ @ @ @ @ @ @   -- conv. layer 5x5x1=>4 stride 1        W1 [5, 5, 1, 4]        B1 [4]
# ∶∶∶∶∶∶∶∶∶∶∶∶∶∶∶∶∶∶∶                                           Y1 [batch, 28, 28, 4]
#   @ @ @ @ @ @ @ @     -- conv. layer 5x5x4=>8 stride 2        W2 [5, 5, 4, 8]        B2 [8]
#   ∶∶∶∶∶∶∶∶∶∶∶∶∶∶∶                                             Y2 [batch, 14, 14, 8]
#     @ @ @ @ @ @       -- conv. layer 4x4x8=>12 stride 2       W3 [4, 4, 8, 12]       B3 [12]
#     ∶∶∶∶∶∶∶∶∶∶∶                                               Y3 [batch, 7, 7, 12] => reshaped to YY [batch, 7*7*12]
#      \x/x\x\x/        -- fully connected layer (relu)         W4 [7*7*12, 200]       B4 [200]
#       · · · ·                                                 Y4 [batch, 200]
#       \x/x\x/         -- fully connected layer (softmax)      W5 [200, 10]           B5 [10]
#        · · ·                                                  Y [batch, 10]

# input X: 28x28 grayscale images, the first dimension (None) will index the images in the mini-batch
X = tf.placeholder(tf.float32, [None, 32, 32, 1])
# correct answers will go here
Y_ = tf.placeholder(tf.float32, [None, 23])
# variable learning rate
lr = tf.placeholder(tf.float32)

# three convolutional layers with their channel counts, and a
# fully connected layer (tha last layer has 10 softmax neurons)
L1 = 64   # 1st convolutional layer output depth
L2 = 64   # 2nd convolutional layer output depth
L3 = 64   # 3rd convolutional layer
L4 = 128  # 4th convolutional layer
L5 = 128  # 5th convolutional layer
L6 = 128  # 6th convolutional layer
L7 = 256  # 7th convolutional layer
L8 = 256  # fully connected layer

W1 = tf.Variable(tf.truncated_normal([3, 3,  1, L1], stddev=0.1))  # 3x3 patch, 1 input channel, L1 output channels
B1 = tf.Variable(tf.ones([L1])/23)

W2 = tf.Variable(tf.truncated_normal([3, 3, L1, L2], stddev=0.1))
B2 = tf.Variable(tf.ones([L2])/23)

W3 = tf.Variable(tf.truncated_normal([3, 3, L2, L3], stddev=0.1))
B3 = tf.Variable(tf.ones([L3])/23)

W4 = tf.Variable(tf.truncated_normal([3, 3, L3, L4], stddev=0.1))
B4 = tf.Variable(tf.ones([L4])/23)

W5 = tf.Variable(tf.truncated_normal([3, 3, L4, L5], stddev=0.1))
B5 = tf.Variable(tf.ones([L5])/23)

W6 = tf.Variable(tf.truncated_normal([3, 3, L5, L6], stddev=0.1))
B6 = tf.Variable(tf.ones([L6])/23)

W7 = tf.Variable(tf.truncated_normal([8 * 8*L6, L7], stddev=0.1))
B7 = tf.Variable(tf.ones([L7])/23)

W8 = tf.Variable(tf.truncated_normal([L7, 23], stddev=0.1))
B8 = tf.Variable(tf.ones([23])/23)

# test flag for batch norm
tst = tf.placeholder(tf.bool)
iter = tf.placeholder(tf.int32)
# dropout probability
pkeep = tf.placeholder(tf.float32)
pkeep_conv = tf.placeholder(tf.float32)

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

# The model
train = tf.placeholder(tf.bool)
pooling_kernel=(1,2,2,1)
Y1c = tf.nn.conv2d(X, W1, strides=[1, 1, 1, 1], padding='SAME')#+B1
Y1n, update_ema1 = batchnorm(Y1c, tst, iter, B1)
Y1r = tf.nn.relu(Y1n)

Y2c = tf.nn.conv2d(Y1r, W2, strides=[1, 1, 1, 1], padding='SAME')#+B2
Y2n, update_ema2 = batchnorm(Y2c, tst, iter, B2)
Y2r = tf.nn.relu(Y2n)

Y3c = tf.nn.conv2d(Y2r, W3, strides=[1, 1, 1, 1], padding='SAME')#+B2
Y3n, update_ema3 = batchnorm(Y3c, tst, iter, B3)
Y3r = tf.nn.relu(Y3n)
Y3p = tf.nn.max_pool(Y3r, ksize=pooling_kernel, strides=pooling_kernel,padding='SAME')

Y4c = tf.nn.conv2d(Y3p, W4, strides=[1, 1, 1, 1], padding='SAME')#+B3
Y4n, update_ema4 = batchnorm(Y4c, tst, iter, B4)
Y4r = tf.nn.relu(Y4n)

Y5c = tf.nn.conv2d(Y4r, W5, strides=[1, 1, 1, 1], padding='SAME')#+B3
Y5n, update_ema5 = batchnorm(Y5c, tst, iter, B5)
Y5r = tf.nn.relu(Y5n)

Y6c = tf.nn.conv2d(Y5r, W6, strides=[1, 1, 1, 1], padding='SAME')#+B3
Y6n, update_ema6 = batchnorm(Y6c, tst, iter, B6)
Y6r = tf.nn.relu(Y6n)
Y6p = tf.nn.max_pool(Y6r, ksize=pooling_kernel, strides=pooling_kernel,padding='SAME')
# reshape the output from the third convolution for the fully connected layer
Y6s = tf.reshape(Y6p, shape=[-1, 8 * 8 * L6])

Y7l = tf.matmul(Y6s, W7)#+B4
Y7n, update_ema7 = batchnorm(Y7l, tst, iter, B7)
Y7r = tf.nn.relu(Y7n)
Ylogits = tf.matmul(Y7r, W8) + B8
Y = tf.nn.softmax(Ylogits)

update_ema = tf.group(update_ema1, update_ema2, update_ema3, update_ema4, update_ema5, update_ema6, update_ema7)

# cross-entropy loss function (= -sum(Y_i * log(Yi)) ), normalised for batches of 100  images
# TensorFlow provides the softmax_cross_entropy_with_logits function to avoid numerical stability
# problems with log(0) which is NaN
cross_entropy = tf.nn.softmax_cross_entropy_with_logits(logits=Ylogits, labels=Y_)
cross_entropy = tf.reduce_mean(cross_entropy)

# accuracy of the trained model, between 0 (worst) and 1 (best)
correct_prediction = tf.equal(tf.argmax(Y, 1), tf.argmax(Y_, 1))
accuracy = tf.reduce_mean(tf.cast(correct_prediction, tf.float32))
predict = Y

# training step, the learning rate is a placeholder
train_step = tf.train.AdamOptimizer(lr).minimize(cross_entropy)

# init
init = tf.global_variables_initializer()
sess = tf.Session()
sess.run(init)

for i in range(500):
  max_learning_rate = 0.003
  min_learning_rate = 0.0001
  decay_speed = 2000.0
  learning_rate = min_learning_rate + (max_learning_rate - min_learning_rate) * math.exp(-i/decay_speed)
  #learning_rate = 0.002
  xs = np.fromfile(fileName+"seisd"+str(i)+".dat", dtype=np.single)
  ys = np.fromfile(fileName+"labeld"+str(i)+".dat", dtype=np.single)
  ys = ys.reshape(100,23)
  xs = xs.reshape(100,32,32)
  xsa = np.ndarray(shape=(100,32,32,1),dtype=np.single)
  n1,n2,n3=32,32,100
  for i3 in range(n3):
    for i2 in range(n2):
      for i1 in range(n1):
        xsa[i3][i2][i1][0] = xs[i3][i2][i1]
  sess.run(train_step, {X: xsa, Y_: ys, lr: learning_rate, tst: False, iter: i})
  sess.run(update_ema, {X: xsa, Y_: ys, tst: False, iter: i})
  if(i%20==0):
    p,a, c = sess.run([predict,accuracy, cross_entropy], {X: txsa, Y_: tys, tst: True})
    p.tofile(fileName+"test"+str(i)+".dat", format="%4")
    print(str(i) + ": accuracy:" + str(a) + " loss: " + str(c))
    k = 0
    #fault = np.ndarray(shape=(102,101),dtype=np.single)
    #dip = np.ndarray(shape=(102,101),dtype=np.single)
    fault = np.ndarray(shape=(1,818),dtype=np.single)
    dip = np.ndarray(shape=(1,818),dtype=np.single)
    for i1 in range(818):
      pd = sess.run(predict, {X: [seista[i1]], tst: True})
      ps = pd[0][0:22]
      im = np.argmax(ps)
      dip[0][i1] = im
      fault[0][i1] = ps[im]
      print max(ps)
      print ps[im]
    fault.tofile(fileName+"faultF3d"+str(i)+".dat", format="%4")
    dip.tofile(fileName+"dipF3d"+str(i)+".dat", format="%4")




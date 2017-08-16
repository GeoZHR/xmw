import tensorflow as tf
import math 

#from tensorflow.examples.tutorials.mnist import input_data
#mnist = input_data.read_data_sets('MNIST_data', one_hot=True, reshape=False)
import numpy as np

fileName = '../../../data/seis/dfault/fake/2d/'
txs = np.fromfile(fileName+"seisd1.dat", dtype=np.single)
tys = np.fromfile(fileName+"labeld1.dat", dtype=np.single)
txs = txs.reshape(100,32,32)
tys = tys.reshape(100,22)
txsa = np.ndarray(shape=(100,32,32,1),dtype=np.single)
n1,n2,n3=32,32,100
for i3 in range(n3):
  for i2 in range(n2):
    for i1 in range(n1):
      txsa[i3][i2][i1][0] = txs[i3][i2][i1]
# input X: 28x28 grayscale images, the first dimension (None) will index the images in the mini-batch
X = tf.placeholder(tf.float32, [None, 32, 32, 1])
# correct answers will go here
Y_ = tf.placeholder(tf.float32, [None, 22])
# variable learning rate
lr = tf.placeholder(tf.float32)
train = tf.placeholder(tf.bool)


# three convolutional layers with their channel counts, and a
# fully connected layer (tha last layer has 10 softmax neurons)
L1 = 64  # first convolutional layer output depth
L2 = 64  # second convolutional layer output depth
L3 = 64  # second convolutional layer output depth
L4 = 128  # second convolutional layer output depth
L5 = 128  # second convolutional layer output depth
L6 = 128  # second convolutional layer output depth
L7 = 256  # third convolutional layer
L8 = 256  # fully connected layer

W1 = tf.Variable(tf.truncated_normal([3, 3,  1, L1], stddev=0.1))  # 3x5 patch, 1 input channel, L1 output channels
B1 = tf.Variable(tf.ones([L1])/22)
W2 = tf.Variable(tf.truncated_normal([3, 3, L1, L2], stddev=0.1))  
B2 = tf.Variable(tf.ones([L2])/22)
W3 = tf.Variable(tf.truncated_normal([3, 3, L2, L3], stddev=0.1))  
B3 = tf.Variable(tf.ones([L3])/22)
W4 = tf.Variable(tf.truncated_normal([3, 3, L3, L4], stddev=0.1))  
B4 = tf.Variable(tf.ones([L4])/22)
W5 = tf.Variable(tf.truncated_normal([3, 3, L4, L5], stddev=0.1))  
B5 = tf.Variable(tf.ones([L5])/22)
W6 = tf.Variable(tf.truncated_normal([3, 3, L5, L6], stddev=0.1))
B6 = tf.Variable(tf.ones([L6])/22)
W7 = tf.Variable(tf.truncated_normal([3, 3, L6, L7], stddev=0.1))
B7 = tf.Variable(tf.ones([L7])/22)
W8 = tf.Variable(tf.truncated_normal([8 * 8*L7, L8], stddev=0.1))
B8 = tf.Variable(tf.ones([L8])/22)
W9 = tf.Variable(tf.truncated_normal([L8, 22], stddev=0.1))
B9 = tf.Variable(tf.ones([22])/22)

# The model
pooling_kernel=(1,2,2,1)
stride = 1  # output is 28x28
Y1  = tf.nn.conv2d( X, W1, strides=[1, stride, stride, 1], padding='SAME') + B1
#Y1n = tf.contrib.layers.batch_norm(Y1, fused=True, is_training=train)
Y1r = tf.nn.relu(Y1)

stride = 1  # outputk is 28x28
Y2 = tf.nn.conv2d(Y1r, W2, strides=[1, stride, stride, 1], padding='SAME') + B2
#Y2n = tf.contrib.layers.batch_norm(Y2, fused=True, is_training=train)
Y2r = tf.nn.relu(Y2)

stride = 1  # output is 28x28
Y3 = tf.nn.conv2d(Y2r, W3, strides=[1, stride, stride, 1], padding='SAME') + B3
#Y3n = tf.contrib.layers.batch_norm(Y3, fused=True, is_training=train)
Y3r = tf.nn.relu(Y3)
Y3p = tf.nn.max_pool(Y3r, ksize=pooling_kernel, strides=pooling_kernel,padding='SAME')

stride = 1  # output is 28x28
Y4 = tf.nn.conv2d(Y3p, W4, strides=[1, stride, stride, 1], padding='SAME') + B4
#Y4n = tf.contrib.layers.batch_norm(Y4, fused=True, is_training=train)
Y4r = tf.nn.relu(Y4)

stride = 1  # output is 28x28
Y5 = tf.nn.conv2d(Y4r, W5, strides=[1, stride, stride, 1], padding='SAME') + B5
#Y5n = tf.contrib.layers.batch_norm(Y5, fused=True, is_training=train)
Y5r = tf.nn.relu(Y5)

stride = 1  # output is 14x14
Y6 = tf.nn.conv2d(Y5r, W6, strides=[1, stride, stride, 1], padding='SAME') + B6
#Y6n = tf.contrib.layers.batch_norm(Y6, fused=True, is_training=train)
Y6r = tf.nn.relu(Y6)

stride = 1  # output is 7x7
Y7 = tf.nn.conv2d(Y6r, W7, strides=[1, stride, stride, 1], padding='SAME') + B7
#Y7n = tf.contrib.layers.batch_norm(Y7, fused=True, is_training=train)
Y7r = tf.nn.relu(Y7)
Y7p = tf.nn.max_pool(Y7r, ksize=pooling_kernel, strides=pooling_kernel,padding='SAME')
# reshape the output from the third convolution for the fully connected layer
YY = tf.reshape(Y7p, shape=[-1, 8*8*L7])

Y8 = tf.nn.relu(tf.matmul(YY, W8) + B8)
Ylogits = tf.matmul(Y8, W9) + B9
Y = tf.nn.softmax(Ylogits)

# cross-entropy loss function (= -sum(Y_i * log(Yi)) ), normalised for batches of 100  images
# TensorFlow provides the softmax_cross_entropy_with_logits function to avoid numerical stability
# problems with log(0) which is NaN
cross_entropy = tf.nn.softmax_cross_entropy_with_logits(logits=Ylogits, labels=Y_)
cross_entropy = tf.reduce_mean(cross_entropy)*10000

# accuracy of the trained model, between 0 (worst) and 1 (best)
predict = Y
correct_prediction = tf.equal(tf.argmax(Y, 1), tf.argmax(Y_, 1))
#correct_prediction = tf.equal(Y, Y_)
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
  decay_speed = 2000.0
  learning_rate = min_learning_rate + (max_learning_rate - min_learning_rate) * math.exp(-i/decay_speed)
  #learning_rate = 0.002
  xs = np.fromfile(fileName+"seisd"+str(i)+".dat", dtype=np.single)
  ys = np.fromfile(fileName+"labeld"+str(i)+".dat", dtype=np.single)
  ys = ys.reshape(100,22)
  xs = xs.reshape(100,32,32)
  xsa = np.ndarray(shape=(100,32,32,1),dtype=np.single)
  n1,n2,n3=32,32,100
  for i3 in range(n3):
    for i2 in range(n2):
      for i1 in range(n1):
        xsa[i3][i2][i1][0] = xs[i3][i2][i1]
  sess.run(train_step, {X: xsa, Y_: ys, lr: learning_rate, train: True})
  if(i%20==0):
    p,a, c = sess.run([predict,accuracy, cross_entropy], {X: xsa, Y_: ys, train: False})
    print p
    p.tofile(fileName+"test"+str(i)+".dat", format="%4")
    print(str(i) + ": accuracy:" + str(a) + " loss: " + str(c))


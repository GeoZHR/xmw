import tensorflow as tf
import math 

#from tensorflow.examples.tutorials.mnist import input_data
#mnist = input_data.read_data_sets('MNIST_data', one_hot=True, reshape=False)
import numpy as np

fileName = '../../../data/seis/dfault/fake/2d/'
txs = np.fromfile(fileName+"attr1.dat", dtype=np.single)
tys = np.fromfile(fileName+"mark1.dat", dtype=np.single)
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


# three convolutional layers with their channel counts, and a
# fully connected layer (tha last layer has 10 softmax neurons)
K = 4  # first convolutional layer output depth
L = 8  # second convolutional layer output depth
M = 12  # third convolutional layer
N = 200  # fully connected layer

W1 = tf.Variable(tf.truncated_normal([5, 5, 1, K], stddev=0.1))  # 5x5 patch, 1 input channel, K output channels
B1 = tf.Variable(tf.ones([K])/22)
W2 = tf.Variable(tf.truncated_normal([5, 5, K, L], stddev=0.1))
B2 = tf.Variable(tf.ones([L])/22)
W3 = tf.Variable(tf.truncated_normal([4, 4, L, M], stddev=0.1))
B3 = tf.Variable(tf.ones([M])/22)

W4 = tf.Variable(tf.truncated_normal([7 * 7 * M, N], stddev=0.1))
B4 = tf.Variable(tf.ones([N])/22)
W5 = tf.Variable(tf.truncated_normal([N, 22], stddev=0.1))
B5 = tf.Variable(tf.ones([22])/22)

# The model
stride = 1  # output is 28x28
Y1 = tf.nn.relu(tf.nn.conv2d(X, W1, strides=[1, stride, stride, 1], padding='SAME') + B1)
stride = 2  # output is 14x14
Y2 = tf.nn.relu(tf.nn.conv2d(Y1, W2, strides=[1, stride, stride, 1], padding='SAME') + B2)
stride = 2  # output is 7x7
Y3 = tf.nn.relu(tf.nn.conv2d(Y2, W3, strides=[1, stride, stride, 1], padding='SAME') + B3)

# reshape the output from the third convolution for the fully connected layer
YY = tf.reshape(Y3, shape=[-1, 7 * 7 * M])

Y4 = tf.nn.relu(tf.matmul(YY, W4) + B4)
Ylogits = tf.matmul(Y4, W5) + B5
Y = tf.nn.softmax(Ylogits)

# cross-entropy loss function (= -sum(Y_i * log(Yi)) ), normalised for batches of 100  images
# TensorFlow provides the softmax_cross_entropy_with_logits function to avoid numerical stability
# problems with log(0) which is NaN
cross_entropy = tf.nn.softmax_cross_entropy_with_logits(logits=Ylogits, labels=Y_)
cross_entropy = tf.reduce_mean(cross_entropy)*100

# accuracy of the trained model, between 0 (worst) and 1 (best)
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
  decay_speed = 2000.0
  learning_rate = min_learning_rate + (max_learning_rate - min_learning_rate) * math.exp(-i/decay_speed)
  xs = np.fromfile(fileName+"attr"+str(i)+".dat", dtype=np.single)
  ys = np.fromfile(fileName+"mark"+str(i)+".dat", dtype=np.single)
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
  sess.run(train_step, {X: xsa, Y_: ys, lr: learning_rate})
  if(i%20==0):
    #a, c = sess.run([accuracy, cross_entropy], {X: batch_xs, Y_: batch_ys})
    a, c = sess.run([accuracy, cross_entropy], {X: xsa, Y_: ys})
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

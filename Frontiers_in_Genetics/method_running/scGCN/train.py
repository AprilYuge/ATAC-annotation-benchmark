from memory_profiler import memory_usage
save = True

import time
import os
import sys
from os import path
from pathlib import Path
import time
import numpy as np
import pickle as pkl
# import tensorflow as tf
import tensorflow.compat.v1 as tf
tf.disable_eager_execution()
from utils import *
from tensorflow.python.saved_model import tag_constants
from models import scGCN
#sys.stdout = open("output_log.txt", "w")

import warnings
warnings.filterwarnings("ignore")
#' del_all_flags(FLAGS)

# Set random seed
seed = 123
np.random.seed(seed)
tf.compat.v1.set_random_seed(seed)
# tf.set_random_seed(seed)

# Get parent folder
import sys
assert len(sys.argv) == 3, "parameters needed: output path, tmp folder"
output_folder = str(sys.argv[1])
result_folder = str(sys.argv[1])
scGCN_folder = str(sys.argv[2])

def main():
    
    start =  time.time()

    # scGCN_folder = path.dirname(path.abspath(__file__))
    Path(output_folder).mkdir(parents=True, exist_ok=True)
    # Settings
    flags = tf.app.flags
    FLAGS = flags.FLAGS
    flags.DEFINE_string('dataset', '%s/input' % scGCN_folder, 'data dir')
    flags.DEFINE_string('output', '%s' % output_folder, 'predicted results')
    flags.DEFINE_bool('graph', True, 'select the optional graph.')
    flags.DEFINE_string('model', 'scGCN','Model string.')
    flags.DEFINE_float('learning_rate', 0.01, 'Initial learning rate.')
    flags.DEFINE_integer('epochs', 200, 'Number of epochs to train.')
    flags.DEFINE_integer('hidden1', 32, 'Number of units in hidden layer 1.')
    #flags.DEFINE_integer('hidden2', 32, 'Number of units in hidden layer 2.')
    flags.DEFINE_float('dropout', 0, 'Dropout rate (1 - keep probability).')
    flags.DEFINE_float('weight_decay', 0,
                       'Weight for L2 loss on embedding matrix.')
    flags.DEFINE_integer('early_stopping', 10,
                         'Tolerance for early stopping (# of epochs).')
    flags.DEFINE_integer('max_degree', 3, 'Maximum Chebyshev polynomial degree.')

    # Load data
    # TODO this `test_on_train` is just a temperary setting, should make explicit interface for train and test dataset input!
    adj, features, labels_binary_train, labels_binary_val, train_mask, pred_mask, val_mask, new_label, true_label, index_guide, rename, data2_index = load_data(
        FLAGS.dataset,rgraph=FLAGS.graph)
    support = [preprocess_adj(adj)]
    num_supports = 1
    model_func = scGCN

    # Define placeholders
    placeholders = {
        'support':
        [tf.sparse_placeholder(tf.float32) for _ in range(num_supports)],
        'features':
        tf.sparse_placeholder(tf.float32,
                              shape=tf.constant(features[2], dtype=tf.int64)),
        'labels':
        tf.placeholder(tf.float32, shape=(None, labels_binary_train.shape[1])),
        'labels_mask':
        tf.placeholder(tf.int32),
        'dropout':
        tf.placeholder_with_default(0., shape=()),
        'num_features_nonzero':
        tf.placeholder(tf.int32)  # helper variable for sparse dropout
    }

    # Create model
    model = model_func(placeholders, input_dim=features[2][1], logging=True)


    # Define model evaluation function
    def evaluate(features, support, labels, mask, placeholders):
        t_test = time.time()
        feed_dict_val = construct_feed_dict(features, support, labels, mask,
                                            placeholders)
        outs_val = sess.run([model.loss, model.accuracy], feed_dict=feed_dict_val)
        return outs_val[0], outs_val[1], (time.time() - t_test)


    # Initialize session
    sess = tf.Session()
    # Init variables
    sess.run(tf.global_variables_initializer())

    train_accuracy = []
    train_loss = []
    val_accuracy = []
    val_loss = []
    # test_accuracy = []
    # test_loss = []

    # Train model

    #configurate checkpoint directory to save intermediate model training weights
    saver = tf.train.Saver()
    save_dir = os.path.join(scGCN_folder, 'checkpoints/')
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    save_path = os.path.join(save_dir, 'best_validation')

    for epoch in range(FLAGS.epochs):
        t = time.time()
        # Construct feed dictionary
        feed_dict = construct_feed_dict(features, support, labels_binary_train,
                                        train_mask, placeholders)
        feed_dict.update({placeholders['dropout']: FLAGS.dropout})
        # Training step
        outs = sess.run([model.opt_op, model.loss, model.accuracy],
                        feed_dict=feed_dict)
        train_accuracy.append(outs[2])
        train_loss.append(outs[1])
        # Validation
        cost, acc, duration = evaluate(features, support, labels_binary_val,
                                       val_mask, placeholders)
        val_loss.append(cost)
        val_accuracy.append(acc)
    #     test_cost, test_acc, test_duration = evaluate(features, support,
    #                                                   labels_binary_test,
    #                                                   test_mask, placeholders)
    #     test_accuracy.append(test_acc)
    #     test_loss.append(test_cost)
        saver.save(sess=sess, save_path=save_path)
        print("Epoch:", '%04d' % (epoch + 1), "train_loss=",
              "{:.5f}".format(outs[1]), "train_acc=", "{:.5f}".format(outs[2]),
              "val_loss=", "{:.5f}".format(cost), "val_acc=", "{:.5f}".format(acc),
              "time=", "{:.5f}".format(time.time() - t))
        if epoch > FLAGS.early_stopping and val_loss[-1] > np.mean(
                val_loss[-(FLAGS.early_stopping + 1):-1]):
            print("Early stopping...")
            break

    print("Finished Training....")

    #'  outputs
    all_mask = np.array([True] * len(train_mask))
    labels_binary_all = new_label

    feed_dict_all = construct_feed_dict(features, support, labels_binary_all,
                                        all_mask, placeholders)
    feed_dict_all.update({placeholders['dropout']: FLAGS.dropout})

    activation_output = sess.run(model.activations, feed_dict=feed_dict_all)[1]
    predict_output = sess.run(model.outputs, feed_dict=feed_dict_all)

    #' accuracy on all masks
    ab = sess.run(tf.nn.softmax(predict_output))
    all_prediction = sess.run(
        tf.equal(sess.run(tf.argmax(ab, 1)),
                 sess.run(tf.argmax(labels_binary_all, 1))))

    #' accuracy on prediction masks
    acc_train = np.sum(all_prediction[train_mask]) / np.sum(train_mask)
    # acc_test = np.sum(all_prediction[test_mask]) / np.sum(test_mask)
    acc_val = np.sum(all_prediction[val_mask]) / np.sum(val_mask)
    # acc_pred = np.sum(all_prediction[pred_mask]) / np.sum(pred_mask)
    print('Checking train/val set accuracy: {}, {}'.format(
        acc_train, acc_val))
    # print('Checking pred set accuracy: {}'.format(acc_pred))

    # save results
    if not os.path.exists(FLAGS.output):
        os.mkdir(FLAGS.output)
    
    if save:    
        # prob matrix
        prob = pd.DataFrame(ab[pred_mask,], index=data2_index, columns=rename.keys())
        prob.to_csv(FLAGS.output+'/prob.csv')
        # predicted label
        pred = prob.idxmax(axis=1)
        pred.to_csv(FLAGS.output+'/pred.csv')

    end = time.time()
    print('Running time: %.2f sec' % (end-start))
    with open(path.join(result_folder, 'time.txt'), "w") as file:
        file.write(str(end-start))
    
peak_mem_usage = memory_usage(main, max_iterations=1, max_usage=True)
print('Peak memory usage: %.2f MB' % peak_mem_usage)
with open(path.join(result_folder, 'memory.txt'), "w") as file:
    file.write(str(peak_mem_usage))

# #' save the predicted labels of query data
# #os.mkdir(FLAGS.output)
# scGCN_all_labels = true_label.values.flatten()  #' ground truth
# np.savetxt(FLAGS.output + '/scGCN_all_input_labels.csv',scGCN_all_labels,delimiter=',',comments='',fmt='%s')
# np.savetxt(FLAGS.output+'/scGCN_query_mask.csv',pred_mask,delimiter=',',comments='',fmt='%s')
# ab = sess.run(tf.nn.softmax(predict_output))
# all_binary_prediction = sess.run(tf.argmax(ab, 1))  #' predict catogrized labels
# all_binary_labels = sess.run(tf.argmax(labels_binary_all, 1))  #' true catogrized labels
# np.savetxt(FLAGS.output+'/scGCN_all_predicted_prob.csv',ab,delimiter=',',comments='',fmt='%f') # yuge
# np.savetxt(FLAGS.output+'/scGCN_all_binary_predicted_labels.csv',all_binary_prediction,delimiter=',',comments='',fmt='%f')
# np.savetxt(FLAGS.output+'/scGCN_index_guide.csv',index_guide,delimiter=',',comments='',fmt='%f')
# np.savetxt(FLAGS.output+'/scGCN_all_binary_input_labels.csv',all_binary_labels,delimiter=',',comments='',fmt='%f')


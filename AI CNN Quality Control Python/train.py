###############################################################
# Treinamento de rede neural convolucional
# Testado com sucesso usando Tensorflow 2.4
###############################################################

import cv2

import tensorflow as tf
from tensorflow.keras import layers
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Conv2D, Flatten, Dropout, MaxPooling2D
from tensorflow.keras.preprocessing.image import ImageDataGenerator
import tensorflow_hub as hub
import PIL.Image as Image
from keras.preprocessing.image import img_to_array
from keras.preprocessing.image import array_to_img
import os
import os.path
from os import path
from pathlib import Path
import pandas as pd
import sys
import fnmatch
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.patches as patches
import numpy as np
from numpy import asarray
import glob

###############################################################


root_dir = "POS"
root_file = root_dir + '/folders.txt'
fid_root = open(root_file, "r")

ok = []
#ko = []
left_corner = []
top_corner = []
right_corner = []
bottom_corner = []
root = []

cnt = 0
while True:
  root_dir = fid_root.readline()
  #print(root_dir)
    
  if root_dir == '':
    break
    
  root_dir = root_dir[0:-1]
  root.append(root_dir)

  ok.append(fid_root.readline()[0:-1])
# ko.append(fid_root.readline()[0:-1])

  print('\nPre-analysis: ' + root_dir)
  print('OK: ', ok[cnt])
#  print('KO: ', ko[cnt])

  crop_file = root_dir + '/crop.txt'
  print('Crop: ', crop_file)
  fid_crop = open(crop_file, "r")
    
  left_corner.append(fid_crop.readline()[0:-1])
  top_corner.append(fid_crop.readline()[0:-1])
  right_corner.append(fid_crop.readline()[0:-1])
  bottom_corner.append(fid_crop.readline()[0:-1])
        
  fid_crop.close()
  
  cnt = cnt + 1

fid_root.close()

print('\ncnt: ', cnt)
print('root: ', root)
print('OK: ', ok)
#print('KO: ', ko)

print('Left: ', left_corner)
print('Top: ', top_corner)
print('Right: ', right_corner)
print('Bottom: ', bottom_corner)

print('root_dir: ', root)

###############################################################

epochs = 3 # 2 # 5 # 10 # 20
steps_per_epoch = 1 # Nao utilizado (automatico)

###############################################################

# Dimensoes fixas para Mobilenet v2 (224 x 224 x 3)
IMG_WIDTH = 224 # 224
IMG_HEIGHT = 224 # 224

# Parametros do treinamento
num_epochs = epochs
steps_per_epoch = steps_per_epoch

# Path das imagens para treinamento/validacao
data_root = 'POS/C680'
data_pre = data_root + '/.pre'
print('data_pre: ' + str(data_pre))

for i in range(0,-1): # !!! Nao executa o loop !!!
#for i in range(0,cnt): # Executa o loop

  test_name = root[i] 
  print('Testing: ' + test_name)

  ### Recorte da imagem
  LEFT_CORNER = left_corner[i]
  TOP_CORNER = top_corner[i]
  RIGHT_CORNER = right_corner[i]
  BOTTOM_CORNER = bottom_corner[i]

###############################################################

  if path.exists(data_pre) == False:
    os.mkdir(data_pre)

  print(data_root)

  final_folder = Path(test_name).parts
  final_folder = final_folder[-1]
  print('final_folder: ' + final_folder)

  if path.exists(data_pre + '/' + final_folder) == False:
    os.mkdir(data_pre + '/' + final_folder)

# !!! Faz o crop e gera a pasta .pre

  files = glob.glob(test_name + '/*.jpg')
  for file in files:
    #print(file)
    image = cv2.imread(file)   
    image = image[int(TOP_CORNER):int(BOTTOM_CORNER), int(LEFT_CORNER):int(RIGHT_CORNER)]
    image = cv2.resize(image, (IMG_WIDTH, IMG_HEIGHT), interpolation = cv2.INTER_AREA)

    final_file = Path(file).parts
    final_file = final_file[-1]   
    #print(data_pre + '/' + final_folder + '/' final_file)
    cv2.imwrite(data_pre + '/' + final_folder + '/'+ final_file, image)

###############################################################

# Normaliza pixels das imagens e le as imagens
# O preprocessing padrao estava colocando os pixels das imagens fora da faixa de valores

###
# Data Augmentation
datagen = ImageDataGenerator(rescale=1./255., validation_split=0.05, width_shift_range=0.1, height_shift_range=0.1, rotation_range=1, brightness_range=[0.9,1.1]) # , fill_mode="wrap")
#datagen = ImageDataGenerator(rescale=1./255., validation_split=0.2)

# !!! Verificar se ainda deve manter class_mode="categorical"

# Path das imagens para treinamento/validacao
data_root = 'POS/C680'
data_pre = data_root + '/.pre'
print('data_pre: ' + str(data_pre))

train_generator = datagen.flow_from_directory( str(data_pre), 
                                         target_size=(IMG_HEIGHT, IMG_WIDTH), 
                                         subset="training",
                                         class_mode="categorical" )
        
valid_generator = datagen.flow_from_directory( str(data_pre), 
                                         target_size=(IMG_HEIGHT, IMG_WIDTH), 
                                         subset="validation",
                                         class_mode="categorical" )

# Lista dimensoes das imagens
for image_train_batch, label_train_batch in train_generator:
  #print("Image batch shape: ", image_train_batch.shape)
  #print("Label batch shape: ", label_train_batch.shape)
  break

# Lista dimensoes das imagens
for image_val_batch, label_val_batch in valid_generator:
  #print("Image batch shape: ", image_val_batch.shape)
  #print("Label batch shape: ", label_val_batch.shape)
  break

# Caminho da rede neural convolucional MobileNet (headless) na Internet (Com KerasLayer)
feature_extractor_url = "https://tfhub.dev/google/tf2-preview/mobilenet_v2/feature_vector/4" #@param {type:"string"}
#feature_extractor_url = "https://tfhub.dev/google/tf2-preview/mobilenet_v2/feature_vector/2" #@param {type:"string"}

feature_extractor_layer = hub.KerasLayer(feature_extractor_url,input_shape=(IMG_HEIGHT,IMG_WIDTH,3)) # Com KerasLayer
#feature_extractor_layer = tf.keras.applications.MobileNetV2(input_shape=(IMG_HEIGHT,IMG_WIDTH,3), classifier_activation=None) # Sem KerasLayer

# Define a partrubber_come da rede que nao podera ser treinada (Feature Extractor)
feature_extractor_layer.trainable = False

print('train_generator.num_classes : ', train_generator.num_classes)

# Cria modelo adicionando nova camada (rede neural convolucional)
model = tf.keras.Sequential([
  feature_extractor_layer,
  layers.Dense(train_generator.num_classes) # !!!
])

#print('Summary')
model.summary()

print('Compile')
# Parametros de otimizacao de modelo (minimizacao do erro)
model.compile(
  optimizer=tf.keras.optimizers.Adam(),
  loss=tf.keras.losses.CategoricalCrossentropy(from_logits=True),
  metrics=['acc'])

# Treina o modelo (rede neural convolucional)
#print(ok)
#print(ko)
#class_weight = {0: ok/(ok+ko), 1: ko/(ok+ko)} # !!!

###############################################################

history = model.fit(train_generator,
                    epochs=num_epochs,
                    validation_data=valid_generator) # !!!

#history = model.fit(train_generator,
#                    epochs=num_epochs,
#                    validation_data=valid_generator,
#                    class_weight=class_weight) # !!!

# Exibe os nomes das classes
class_names = sorted(train_generator.class_indices.items(), key=lambda pair:pair[1])
class_names = np.array([key.title() for key, value in class_names])
print('class_names: ',  class_names)

###############################################################

# Classifica imagens do dataset (train)
predicted_batch = model.predict(image_train_batch)
predicted_id = np.argmax(predicted_batch, axis=-1)
predicted_label_batch = class_names[predicted_id]
label_id = np.argmax(label_train_batch, axis=-1)

# Exibe painel com imagens
plt.figure(figsize=(10,9))
plt.subplots_adjust(hspace=0.5)
for n in range(25):
  plt.subplot(5,5,n+1)
  plt.imshow(image_train_batch[n])
  color = "green" if predicted_id[n] == label_id[n] else "red"
  plt.title(predicted_label_batch[n].title(), color=color)
  plt.axis('off')
_ = plt.suptitle("Model predictions (green: correct, red: incorrect)")

plt.show()
#plt.savefig('all_train.png')
#plt.savefig(class_file + '_train.png')
plt.close()

###############################################################

# Classifica imagens do dataset (validation)
predicted_batch = model.predict(image_val_batch)
predicted_id = np.argmax(predicted_batch, axis=-1)
predicted_label_batch = class_names[predicted_id]
label_id = np.argmax(label_val_batch, axis=-1)

# Exibe painel com imagens
plt.figure(figsize=(10,9))
plt.subplots_adjust(hspace=0.5)
for n in range(25):
  plt.subplot(5,5,n+1)
  plt.imshow(image_val_batch[n])
  color = "green" if predicted_id[n] == label_id[n] else "red"
  plt.title(predicted_label_batch[n].title(), color=color)
  plt.axis('off')
_ = plt.suptitle("Model validations (green: correct, red: incorrect)")

plt.show()
#plt.savefig('all_val.png')
#plt.savefig(class_file + '_val.png')
plt.close()

###############################################################

### Salva modelo
test_name = 'POS/C680/all'
export_path = test_name
model.save(export_path + '.h5')

# Salva pasta com o modelo
model = tf.keras.models.load_model(export_path + '.h5', custom_objects={'KerasLayer':hub.KerasLayer}) # Com KerasLayer
tf.saved_model.save(model, export_path + '_frozen')


import os

#from tensorflow.python.compiler.tensorrt import trt_convert as trt
import tensorflow as tf
from tensorflow.keras import layers
import tensorflow_hub as hub
from tensorflow.keras.models import load_model
from tensorflow.keras.preprocessing import image

import keras
from keras.models import load_model

import PIL.Image as Image
from pathlib import Path
from matplotlib import pyplot as plt
import numpy as np
import argparse
import math
import time
import glob
import cv2
import sys

# Classes da rede neural
class_names = ['card_ko', 'card_ok', 'conector_bateria_ko', 'conector_bateria_ok', 'flat_ko', 'flat_ok', 'inserto_ko', 'inserto_ok', 'j15_ko', 'j15_ok', 'j1_ko', 'j1_ok', 'rebarba_ko', 'rebarba_ok', 'rubber_ko', 'rubber_ok']
print(class_names)

# Caminho do diretorio onde sera feita a analise
path_base = 'teste'

# Show graphics?
GRAPHICS = 1
#GRAPHICS = 0

# Dimensoes fixas para Mobilenet v2 (224 x 224 x 3)
IMG_WIDTH = 224 # 224
IMG_HEIGHT = 224 # 224

# Tamanho original da imagem
ORG_WIDTH = 2590
ORG_HEIGHT = 1942

root_dir = 'POS'
POS_model = 'C680'
root_folder = root_dir + '/' + POS_model
print(root_folder)
model_path = root_folder + '/all_frozen'
print(model_path)

folders = glob.glob(root_folder + '/*_ok/')
print(folders)

#############################################################################
# POS
#   models
#############################################################################

files = [f for f in os.listdir(path_base) if (os.path.isfile(os.path.join(path_base,f)))]
#print(files)
lf = len(files)
#print(lf)

#############################################################################

start_time = time.time()

print('Model path: ' + model_path)

# Load converted model
model = tf.keras.models.load_model(model_path)

print("--- Neural Networks Loading takes %s seconds ---" % (time.time() - start_time))

#############################################################################

for ff in range(0,lf):
    
    file = os.path.join(path_base, files[ff])
    print(file)

    print('\nBoard image: ' + file + '\n')

    img_org = cv2.imread(file)
    height, width = img_org.shape[:2]
    #print(width)
    #print(height)

    fac = 3
    w = int(width/fac)
    h = int(height/fac)
    d = (w, h)
    img_show = cv2.resize(img_org, d)

    if (GRAPHICS == 1):
        cv2.imshow(file, img_show)
        #cv2.namedWindow(file, cv2.WINDOW_AUTOSIZE)
        cv2.waitKey(100)
        #cv2.destroyAllWindows()

    #################################################

    start_time = time.time()

    cnt = 0
    for folder in folders:

        #print(folder)

        crop_file = folder + 'crop.txt'
        #print(crop_file)
        fid_crop = open(crop_file, "r")

        LEFT_CORNER = int(fid_crop.readline()[0:-1])
        TOP_CORNER = int(fid_crop.readline()[0:-1])
        RIGHT_CORNER = int(fid_crop.readline()[0:-1])
        BOTTOM_CORNER = int(fid_crop.readline()[0:-1])
        print(LEFT_CORNER)
        print(TOP_CORNER)
        print(RIGHT_CORNER)
        print(BOTTOM_CORNER)

        fid_crop.close()

        h5_file = folder[0:-1] + '.h5'
        #pb_file = folder[0:-1] + '.pb'
        #print(h5_file)

        #################################################

        (head, tail) = os.path.split(folder[0:-1])
        print('*Testing: ' + tail)
        print(file)

        #################################################

        cv_img = cv2.imread(file)
        cv_img = cv_img[int(TOP_CORNER):int(BOTTOM_CORNER), int(LEFT_CORNER):int(RIGHT_CORNER)]
        cv_img = cv2.resize(cv_img, (IMG_WIDTH, IMG_HEIGHT), interpolation = cv2.INTER_AREA)
        final_file = Path(file).parts
        final_file = final_file[-1]
        #print('crops/' + tail + '_' + final_file)
        cv2.imwrite('crops/' + tail + '_' + final_file, cv_img)

        #################################################

        img = image.load_img(file)
        img = img.crop( (LEFT_CORNER, TOP_CORNER, RIGHT_CORNER, BOTTOM_CORNER) )
        img = img.resize( (IMG_WIDTH, IMG_HEIGHT), Image.ANTIALIAS)
        #img.show()
        img_array = image.img_to_array(img)
        img_array = np.expand_dims(img_array, axis=0)

        #################################################

        print('***Model: ')
        print(model)

        current_time = time.time()

        output = model.predict(img_array / 255.)
        print(output)

        lc = len(class_names)
        print(tail)
        for i in range(0,lc):
            if class_names[i] == tail:
                print(class_names[i])
                out = output[0][i]
                print(out)      
 
        print ('Time to predict a single model: %s seconds' % (time.time()-current_time))
        
        #################################################

        cnt = cnt + 1

        if out > 0:
        #if out<0.5:
            color = (0,255,0)
            print("OK!")
        else:
            color = (0,0,255)
            print("KO!")

        print(' ')

        if (GRAPHICS == 1):
            start_point = (int(LEFT_CORNER/fac), int(TOP_CORNER/fac))
            end_point = (int(RIGHT_CORNER/fac), int(BOTTOM_CORNER/fac))
            thickness = 5
            cv2.rectangle(img_show, start_point, end_point, color, thickness)
            cv2.putText(img_show, tail, (int(LEFT_CORNER/fac), int(TOP_CORNER/fac)), cv2.FONT_HERSHEY_SIMPLEX, 0.3, (255, 255, 255), 2, cv2.LINE_AA)
            cv2.putText(img_show, tail, (int(LEFT_CORNER/fac), int(TOP_CORNER/fac)), cv2.FONT_HERSHEY_SIMPLEX, 0.3, (0, 0, 0), 1, cv2.LINE_AA)
            cv2.imshow(file, img_show)
            cv2.waitKey(100)

    print("--- Neural Networks Prediction takes %s seconds ---" % (time.time() - start_time))

    if (GRAPHICS == 1):
        cv2.waitKey(0)

import cv2 as cv
import numpy as np
img= cv.imread('cat1.jpg')
def translate(img,x,y):
    transMat=np.float32([[0,1,x],[0,1,y]])
    dimensions=(img.shape[1],img.shape[0])
    return cv.warpAffine(img,transMat,dimensions)
#-x-left,x-right
#-y-up, y-down
translated=translate(img,100,100)
cv.imshow("Translated",translated)
#rotation
def rotation(img,angle,rotPoint=None):
    (height,width)=img.shape[:2]
    if rotPoint is None:
        rotPoint=(width//2,height//2)
    rotMat=cv.getRotationMatrix2D(rotPoint,angle,1.0)
    dimensions=(width,height)
    return cv.warpAffine(img,rotMat,dimensions)
rotated=rotation(img,45)
cv.imshow("Rotated",rotated)
cv.waitKey(0)



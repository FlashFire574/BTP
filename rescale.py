import cv2 as cv
img= cv.imread('cat1.jpg')
cv.imshow('Cat1',img)

def rescaledimage(img, scale=0.75):
    w = int(img.shape[1]*scale)
    h = int(img.shape[0]*scale)
    dimensions=(w,h)
    return cv.resize(img, dimensions,interpolation=cv.INTER_AREA)

new_img=rescaledimage(img)
cv.imshow('Cat1 small',new_img)
cv.waitKey(0)
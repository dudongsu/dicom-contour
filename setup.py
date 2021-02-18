from setuptools import setup

setup(
     name='dicom-contour',
     packages=['dicom_contour'],    
     version=2.3,
     description='A library which converts DICOM images and contours into numpy arrays. An automated way of extracting image and mask voxels.',
     author='Dongsu Du',
     author_email='dudongsu@gmail.com',
     license='MIT',
     url='https://github.com/dudongsu/dicom-contour',
     install_requires=
	['pydicom', 'numpy', 'matplotlib'],
     keywords=['dicom', 'contour', 'mask', 'medical image'],
     
)

cmake_minimum_required(VERSION 2.4)

PROJECT( GaussianInterpolation )

#The following lines are required to use Dart
ENABLE_TESTING()
INCLUDE(Dart)

# Set up ITK
IF(USE_ITK)
  FIND_PACKAGE(ITK)
  IF(ITK_FOUND)
    INCLUDE(${ITK_USE_FILE})
  ELSE(ITK_FOUND)
    MESSAGE(FATAL_ERROR
      "Cannot build without ITK.  Please set ITK_DIR.")
  ENDIF(ITK_FOUND)
ENDIF(USE_ITK)

INCLUDE (${CMAKE_ROOT}/Modules/FindITK.cmake) 
IF (USE_ITK_FILE) 
  INCLUDE(${USE_ITK_FILE}) 
ENDIF(USE_ITK_FILE) 

LINK_DIRECTORIES(${ITK_LIBRARY_PATH})

IF(BUILD_TESTING)
# non-templated class -- this should be stored in a library and linked in...
ADD_EXECUTABLE(itkGaussianInterpolateImageFunctionTest itkGaussianInterpolateImageFunctionTest.cxx)
TARGET_LINK_LIBRARIES(itkGaussianInterpolateImageFunctionTest ITKIO)

ENDIF(BUILD_TESTING)

###
#  Perform testing
###

ADD_TEST( NN_INTERP ${CMAKE_BINARY_DIR}/itkGaussianInterpolateImageFunctionTest 2
          ${CMAKE_SOURCE_DIR}/r16slice.nii.gz
          ${CMAKE_BINARY_DIR}/r16slice_interp_nn.nii.gz
          50x50 1 1 )
ADD_TEST( LINEAR_INTERP ${CMAKE_BINARY_DIR}/itkGaussianInterpolateImageFunctionTest 2
          ${CMAKE_SOURCE_DIR}/r16slice.nii.gz
          ${CMAKE_BINARY_DIR}/r16slice_interp_linear.nii.gz
          50x50 1 0 )
ADD_TEST( GAUSS1x1_INTERP ${CMAKE_BINARY_DIR}/itkGaussianInterpolateImageFunctionTest 2
          ${CMAKE_SOURCE_DIR}/r16slice.nii.gz
          ${CMAKE_BINARY_DIR}/r16slice_interp_guassian1x1.nii.gz
          50x50 1 2 1x1 3.0 )
ADD_TEST( GAUSS5x5_INTERP ${CMAKE_BINARY_DIR}/itkGaussianInterpolateImageFunctionTest 2
          ${CMAKE_SOURCE_DIR}/r16slice.nii.gz
          ${CMAKE_BINARY_DIR}/r16slice_interp_guassian5x5.nii.gz
          50x50 1 2 5x5 3.0 )
ADD_TEST( GAUSS1x5_INTERP ${CMAKE_BINARY_DIR}/itkGaussianInterpolateImageFunctionTest 2
          ${CMAKE_SOURCE_DIR}/r16slice.nii.gz
          ${CMAKE_BINARY_DIR}/r16slice_interp_guassian1x5.nii.gz
          50x50 1 2 1x5 3.0 )

#include <math.h>
#include "../include/Mixer3D.h"
#include "../include/World.h"
#include "../include/fft.h"
#include "../include/mit_hrtf_lib.h"
#include "../include/complextype.h"


Mixer3D::Mixer3D(int bufSize, int smpRate, int bitD, World *w) :
myWorld(w), bufferSize(bufSize), sampleRate(smpRate), bitDepth(bitD), player(myWorld->getPlayer())
{
    inputAO = new Complex[2*bufferSize];
    overlapInput = new Complex[2 * bufferSize];
    fInput = new Complex[2 * bufferSize];
    fFilter = new Complex[2 * bufferSize];
    
    azimuths = new int[World::MAX_OBJ];
    elevations = new int[World::MAX_OBJ];
    prevAzimuths = new int[World::MAX_OBJ];
    prevElevations = new int[World::MAX_OBJ];

    updateAngles(); // initialize azimuths and elevations

    outputLeft = new Complex*[World::MAX_OBJ];
    outputRight = new Complex*[World::MAX_OBJ];
    overlapLeft = new Complex*[World::MAX_OBJ];
    overlapRight = new Complex*[World::MAX_OBJ];

    leftFilter = new short[2 * bufferSize];
    rightFilter = new short[2 * bufferSize];
    ComplexLeftFilter = new Complex*[World::MAX_OBJ];
    ComplexRightFilter = new Complex*[World::MAX_OBJ];

    for (int i = 0; i < World::MAX_OBJ; i++)
    {
        // TODO: Should we initialize filter values?
        overlapLeft[i] = new Complex[bufferSize];
        overlapRight[i] = new Complex[bufferSize];
        outputLeft[i] = new Complex[2*bufferSize];
        outputRight[i] = new Complex[2*bufferSize];
        ComplexLeftFilter[i] = new Complex[2*bufferSize];
        ComplexRightFilter[i] = new Complex[2*bufferSize];
        
        // TODO: Is there a more efficient way to know which objects
        // have changed position relative to the player?
        // Store some state each time an AudioObject is moved?
        // What if the player moves in the same direction?
        prevAzimuths[i] = INVALID_ANGLE;
        prevElevations[i] = INVALID_ANGLE;
    }    
}

Mixer3D::~Mixer3D()
{
    delete [] inputAO;
    delete [] outputLeft;
    delete [] outputLeft;
    delete [] leftFilter;
    delete [] rightFilter;
    delete [] ComplexLeftFilter;
    delete [] ComplexRightFilter;
    delete [] overlapLeft;
    delete [] overlapRight;
    delete [] overlapInput;
    delete [] prevAzimuths;
    delete [] prevElevations;
    delete [] azimuths;
    delete [] elevations;
}

void Mixer3D::updateAngles() {
    AudioObj* iAudioObj;
    for (int i = 0; i < myWorld->getNumObj(); i++) {
        iAudioObj = myWorld->getAudioObj(i);
        azimuths[i] = player.computeAzimuth(iAudioObj);
        elevations[i] = player.computeElevation(iAudioObj);
    }
}

bool Mixer3D::isPowerOfTwo(int x) {
    return !(x == 0) && !(x & (x - 1));
}

int Mixer3D::loadHRTF(int* pAzimuth, int* pElevation, unsigned int samplerate, unsigned int diffused, Complex *&leftFilterIn, Complex *&rightFilterIn)
{
    int size = mit_hrtf_get(pAzimuth, pElevation, samplerate, diffused, leftFilter, rightFilter);

    for (int i = 0; i < size; i++)
    {
        leftFilterIn[i] = (double)(leftFilter[i]);
        rightFilterIn[i] = (double)(rightFilter[i]);
    }

    return size;
}
void Mixer3D::convolution(Complex *input, Complex *filter, Complex *output, long nSig, long nFil, long nFFT) {

	if (input == NULL || filter == NULL) {
		throw invalid_argument("Input and Filter must be non-NULL.");
	}

	if (!isPowerOfTwo(nFFT) || nFFT < nSig || nFFT < nFil) {
        throw invalid_argument("NFFT must be a power of two, bigger than the signal length, and bigger than the filter length.");
	}

	// Perform FFT on both input and filter.
    // TODO: "Streamline" CFFT class?
	CFFT::Forward(input, fInput, (unsigned int)nFFT);
	CFFT::Forward(filter, fFilter, (unsigned int)nFFT);

	for (int i = 0; i < nFFT; i++) {
		output[i] = fInput[i] * fFilter[i];
    }// SIMD?
	CFFT::Inverse(output, (unsigned int)nFFT);
}

void Mixer3D::stereoConvolution(Complex *input, Complex *leftFilter, Complex *rightFilter, Complex *leftOutput, Complex *rightOutput, long nSIG, long nFIL, long nFFT)
{
    // TODO: Modify parameter name, input is a data member
    convolution(input, leftFilter, leftOutput, nSIG, nFIL, nFFT);
    convolution(input, rightFilter, rightOutput, nSIG, nFIL, nFFT);
}

void Mixer3D::performMix(short *ioDataLeft, short *ioDataRight)
{
    //Zero the output arrays.
    for(int i = 0; i < bufferSize; i++)
    {
       ioDataLeft[i] = 0;
       ioDataRight[i] = 0;
    }// SIMD?
    
    //Iterate through all audio objects, obtain input data and calculate resulting audio data for each object.
    // TODO: Possible bug? Cuts short if one audioObj is inactive?
    AudioObj* iAudioObj;
    float iVolume, iDistance, iAmplitudeFactor;
    for (int i = 0; i < myWorld->getNumObj() && myWorld->getAudioObj(i)->isActive(); i++) {
        // dynamically calculate the Azimuth and Elevation between every object and the player
        updateAngles();
        iAudioObj = myWorld->getAudioObj(i);
   
        // loading in input data for the iteration accordingly
        if (!(iAudioObj->fillAudioData(inputAO, bufferSize))) {
            continue;
        }        
       
        // scale volume according to distance of object from player
        // TODO: Is this realistic?
        iVolume = iAudioObj->getVolume();
        iDistance = player.computeDistanceFrom(iAudioObj) ;
        iAmplitudeFactor = iVolume / iDistance;
        //iAmplitudeFactor = (iVolume / iDistance) / myWorld->getNumObj();
        
        for(int j = 0; j < bufferSize * 2; j++) {
            if ( j >= bufferSize ) {
                // zero pad
                inputAO[j] = 0;
            } else {
                inputAO[j] *= iAmplitudeFactor;
            }
        }// SIMD?
       
        if (azimuths[i] != prevAzimuths[i] ||
           elevations[i] != prevElevations[i]) {
            // object location relative to player has changed, so fetch a new filter
            filterLength = loadHRTF(&azimuths[i], &elevations[i], sampleRate, 1, ComplexLeftFilter[i], ComplexRightFilter[i]);
            
            // zero pad
            for (int j = filterLength; j < 2 * bufferSize; j++) {
                ComplexLeftFilter[i][j] = 0;
                ComplexRightFilter[i][j] = 0;
            }// SIMD?
            
            if (azimuths[i] < 0) {
                   azimuths[i] = -azimuths[i];
            }
            
            
            //Since the filter changed, we perform a convolution on the old input with the new filter and pull out its tail.
            
            stereoConvolution(overlapInput, ComplexLeftFilter[i], ComplexRightFilter[i], outputLeft[i], outputRight[i], bufferSize, filterLength, 2 * bufferSize);
        
            // update the overlap part for the next iteration
            for (int j = 0; j < bufferSize; j++) {
                    overlapLeft[i][j] = outputLeft[i][j + bufferSize];
                    overlapRight[i][j] = outputRight[i][j + bufferSize];
            }
        }// SIMD?
        
        //Perform the convolution of the current input and current filter.
        stereoConvolution(inputAO, ComplexLeftFilter[i], ComplexRightFilter[i], outputLeft[i], outputRight[i], bufferSize, filterLength, 2 * bufferSize);
        
        
        for (int j = 0; j < bufferSize; j++) {
            ioDataLeft[j] +=  (short((outputLeft[i][j].real())/2 + (overlapLeft[i][j].real())/2))/myWorld->getNumObj();
            ioDataRight[j] += (short((outputRight[i][j].real())/2 + (overlapRight[i][j].real())/2))/myWorld->getNumObj();
        }// SIMD?
        
  
        //Output the data to the output arrays in short integer format. In addition to pushing the output of the main
        //convolution, we also need to add the overlapped tail of the last output and divide by 2. Finally, we need
        //to divide by the number of audio objects to ensure no clipping.
        for (int j = 0; j < bufferSize; j++)
        {
            
            // this one has logical error
//            ioDataLeft[j]  += (short)( (outputLeft[i][j].real() + overlapLeft[i][j].real()) / 2*myWorld->getNumObj() );
//            ioDataRight[j] += (short)( (outputRight[i][j].real() + overlapRight[i][j].real()) / 2*myWorld->getNumObj() );
            
            // this one does not. 
            ioDataLeft[j]  += (short)( (outputLeft[i][j].real() + overlapLeft[i][j].real()) / (2*myWorld->getNumObj()) );
            ioDataRight[j] += (short)( (outputRight[i][j].real() + overlapRight[i][j].real()) / (2*myWorld->getNumObj()) );
        }// SIMD?
        
        //Updating the overlapInput for the next iteration for the correpsonding obejct
        for (int j = 0; j < bufferSize * 2; j++) {
            if (j >= bufferSize) {
                overlapInput[j] = 0;
            } else {
                overlapInput[j] = inputAO[j];
            }
        }// SIMD?

        // TODO: If the filter has been changed, didn't we already do this?
        //Updating the default overlap information for the next iteration if the filter won't be changed
        for (int j = 0 ; j < bufferSize; j++) {
            overlapLeft[i][j] = outputLeft[i][j + bufferSize];
            overlapRight[i][j] = outputRight[i][j + bufferSize];
        }// SIMD?
    
        //storing the Azimuth value in this iteration for the comparison for the next iteration so that
        //we can know that whether the filter needs to be changed in the next iteration.
        prevAzimuths[i] = azimuths[i];
    }
}

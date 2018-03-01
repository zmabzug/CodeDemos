## Data-Processing-Pipeline

Welcome! This file contains code samples in MATLAB and the [Spike2 scripting language](http://ced.co.uk/products/spkdpsl) that I needed to write in order to integrate new hardware and software technologies into our existing data acquisition and processing pipeline. Unfortunately, all of these files require either specialized software or large numbers of supporting functions in order to run, so they are for viewing only. If you'd like to see a demo of this code in action, please contact me!

### _intanBatchImport.s2s_
Previously, a single stream of analog neurophysiological data was recorded in Spike2, and then processed (spike-sorted) in Spike2. In our new data acquisition pipeline, we recorded multiple streams of data simultaneously using a [different piece of software](http://intantech.com/downloads.html), but still wanted to process data in Spike2. Therefore, I had to write a Spike2 script to take a folder of simultaneously-recorded data streams, convert them to Spike2's desired format, apply digital filtering and DC offset correction, and pair each individual data stream with the "trigger" data (see below for more details). This was my first (and only) piece of code written in the Spike2 scripting language, and I had to write it without any formal training in the language.

### _compareTriggers.m_
In our typical data acquisition pipeline, behavioral and neurophysiological data are recorded in parallel on separate hardware systems. To ensure that these two data streams can be aligned correctly in time, we send "trigger" pulses from the behavioral recording system to the neural recording system at pre-determined times. By aligning when triggers are sent and received, we can reliably align the timing of our two data streams.

However, I discovered that our new neural data acquisition system was prone to receiving extra "phantom" triggers, leading to mismatch between the number of triggers sent and received. I wrote this file to characterize and visualize the relative timing of the sent and received triggers, in order to create a function that would be able to reliably diagnose and remove phantom triggers during data processing.

### _LoadRex\_MergeSpk2\_FixTrigs.m_
This function is responsible for (amongst other things) identifying and removing phantom triggers. It is modified from and built upon code that many people have written over many years, but the relevant portion (which represents some, but not all of my contribution to this function) is contained in lines 60-120.

This code segment extracts the sent and received triggers and steps through them individually, comparing the deviation between each pair. Small deviations (for instance, due to drift) are expected and acceptable, and are (temporarily) adjusted for to prevent large deviations from accumulating. Large deviations, however, are indicative of phantom triggers that need to be removed.

### _rex\_process\_Sp2.m_
Lastly, this function is the only one that a typical user will need to interact with in order to fully process their data into a form usable for analysis. Users can call it from the command line and input their home directory (dependent on the computer), the directory their data is stored in (dependent on the data source), and the files they'd like to process. This file uses _LoadRex\_MergeSpk2\_FixTrigs.m_ to remove any errant triggers, and packages data into matrices for easy analysis. 

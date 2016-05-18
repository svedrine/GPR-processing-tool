# Ascan v.2.0

import os, argparse
import h5py
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from gprMax.exceptions import CmdInputError
from gprMax.receivers import Rx

"""Plots electric and magnetic fields and currents from all receiver points in the given output file. Each receiver point is plotted in a new figure window."""

# Parse command line arguments
parser = argparse.ArgumentParser(description='Plots electric and magnetic fields and currents from all receiver points in the given output file. Each receiver point is plotted in a new figure window.', usage='cd gprMax; python -m tools.plot_Ascan outputfile')
parser.add_argument('outputfile', help='name of output file including path')
parser.add_argument('--outputs', help='outputs to be plotted (Ex, Ey, Ez, Hx, Hy, Hz, Ix, Iy, Iz)', default=Rx.availableoutputs, nargs='+')
parser.add_argument('-fft', action='store_true', default=False, help='plot FFT (single output must be specified)')
parser.add_argument('-gain', choices = ['lin', 'exp'])
parser.add_argument('-p', '-parameters', type= float, help= '(t0, coef)', nargs = 2)
args = parser.parse_args()
	

# Open output file and read some attributes
f = h5py.File(args.outputfile, 'r')
nrx = f.attrs['nrx']
dt = f.attrs['dt']
iterations = f.attrs['Iterations']
time = np.linspace(0, 1, iterations)
time *= (iterations * dt)

# Check there are any receivers
if nrx == 0:
    raise CmdInputError('No receivers found in {}'.format(args.outputfile))

# Check for single output component when doing a FFT
if args.fft:
    if not len(args.outputs) == 1:
        raise CmdInputError('A single output must be specified when using the -fft option')

# New plot for each receiver
for rx in range(1, nrx + 1):
	path = '/rxs/rx%s/' % str(rx)
	availableoutputs = list(f[path].keys())
	
	# If Gain is turned on 
	if args.gain:
		if args.fft:
			raise CmdInputError('Cannot plot FFT with -gain turned on')

		# Create the gained output file on first itteration 
		if rx == 1:
			outputfile = args.outputfile.replace('.out', '_gained.out')
			fg = h5py.File(outputfile, 'w')
		
			# Init datasets' output file
			fg.attrs['Iterations'] = iterations
			fg.attrs['dt'] = dt
			fg.attrs['nrx'] = nrx

			# Find initial time idexation 
			tindex = 0
			tinit = args.p[0]*1e-9
			while time[tindex] <= tinit:
				tindex += 1
			print(tindex)
			# Write the gain function to be applied
			fgain = [1]*len(time[:tindex])
			coef = args.p[1]*1e9
	
			if args.gain == 'lin':
				for t in time[tindex:]:
					fgain.append(coef * (t-tinit) + 1) 
			 
		# Write properties for all availiable outputs
		grp = fg.create_group(path)
		for output in availableoutputs:
			# Initialize
			grp.create_dataset(output, (iterations,), dtype=f['%s%s' %(path, output)].dtype)
			fg['%s%s' %(path, output)][:,] = f['%s%s' %(path, output)][:,]
			# Apply gain
			fg['%s%s' %(path, output)][tindex:,] *= fgain[tindex:]
		
		# Check for single output component when plotting 
		if not len(args.outputs) == 1:
			raise CmdInputError('A single output must be specified when using the --gain option')
		
		# Check for polarity of the single output to be plotted  
		output = args.outputs[0]
		if output[0] == 'm':
			polarity = -1
			outputtext = '-' + output[1:]
			output = output[1:]
		else:
			polarity = 1
			outputtext = output

		# Plot the single output
		data = f['%s%s' %(path, output)][:,] * polarity
		gdata = fg['%s%s' %(path, output)][:,] * polarity
		
		fig, ax = plt.subplots(subplot_kw=dict(xlabel='Time [s]'), num='rx' + str(rx), figsize=(20, 10), facecolor='w', edgecolor='w')
		gs = gridspec.GridSpec(3, 1, hspace=0.3, wspace=0.3)
		
		ax = plt.subplot(gs[0, 0])
		line = ax.plot(time, data,'r', lw=2, label=outputtext)
		ax.set_ylabel('%s, field strength [V/m]' % outputtext)
		
		ax = plt.subplot(gs[1, 0])
		ax.plot(time, fgain,'black', lw=2, label='Fonction gain')
		ax.set_ylabel('Function gain [None]')

		ax = plt.subplot(gs[2, 0])
		gline = ax.plot(time, gdata,'r', lw=2, label=outputtext)
		ax.set_ylabel('%s, field strength [V/m]' % outputtext)
		
		if 'H' in output:
			plt.setp(line, color='g')
			plt.setp(gline, color='g')
			plt.setp(ax, ylabel='%s, field strength [A/m]' % outputtext)
		elif 'I' in output:
			plt.setp(line, color='b')
			plt.setp(gline, color='b')
			plt.setp(ax, ylabel='%s, current [A]' % outputtext)
		
		for ax in fig.axes:
				ax.set_xlim([0, np.amax(time)])
				ax.grid()


	# If gain is turned off  
	else: 		
		# If only a single output is required, create one subplot
		if len(args.outputs) == 1:
		    
		    # Check for polarity of output and if requested output is in file
			if args.outputs[0][0] == 'm':
				polarity = -1
				outputtext = '-' + args.outputs[0][1:]
				output = args.outputs[0][1:]
			else:
				polarity = 1
				outputtext = args.outputs[0]
				output = args.outputs[0]
		    
			if output not in availableoutputs:
				raise CmdInputError('{} output requested to plot, but the available output for receiver 1 is {}'.format(output, ', '.join(availableoutputs)))
		    
			outputdata = f['%s%s' %(path, output)][:] * polarity

		    # Plotting if FFT required
			if args.fft:
				# Calculate magnitude of frequency spectra of waveform
				power = 10 * np.log10(np.abs(np.fft.fft(outputdata))**2)
				freqs = np.fft.fftfreq(power.size, d=dt)

		        # Shift powers so that frequency with maximum power is at zero decibels
				power -= np.amax(power)

		        # Set plotting range to -60dB from maximum power
				pltrange = np.where((np.amax(power[1::]) - power[1::]) > 60)[0][0] + 1
		        # To a maximum frequency
		        #pltrange = np.where(freqs > 2e9)[0][0]
				pltrange = np.s_[0:pltrange]

		        # Plot time history of output component
				fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, num='rx' + str(rx), figsize=(20, 10), facecolor='w', edgecolor='w')
				line1 = ax1.plot(time, outputdata, 'r', lw=2, label=outputtext)
				ax1.set_xlabel('Time [s]')
				ax1.set_ylabel('%s field strength [V/m]' % outputtext)
				ax1.set_xlim([0, np.amax(time)])
				ax1.grid()

		        # Plot frequency spectra
				markerline, stemlines, baseline = ax2.stem(freqs[pltrange], power[pltrange], '-.')
				plt.setp(baseline, 'linewidth', 0)
				plt.setp(stemlines, 'color', 'r')
				plt.setp(markerline, 'markerfacecolor', 'r', 'markeredgecolor', 'r')
				line2 = ax2.plot(freqs[pltrange], power[pltrange], 'r', lw=2)
				ax2.set_xlabel('Frequency [Hz]')
				ax2.set_ylabel('Power [dB]')
				ax2.grid()
		        
		        # Change colours and labels for magnetic field components or currents
				if 'H' in args.outputs[0]:
					plt.setp(line1, color='g')
					plt.setp(line2, color='g')
					plt.setp(ax1, ylabel='%s field strength [A/m]')
					plt.setp(stemlines, 'color', 'g')
					plt.setp(markerline, 'markerfacecolor', 'g', 'markeredgecolor', 'g')
				elif 'I' in args.outputs[0]:
					plt.setp(line1, color='b')
					plt.setp(line2, color='b')
					plt.setp(ax1, ylabel='%s current [A]')
					plt.setp(stemlines, 'color', 'b')
					plt.setp(markerline, 'markerfacecolor', 'b', 'markeredgecolor', 'b')
		        
				plt.show()
		
		    # Plotting if no FFT required
			else:
				fig, ax = plt.subplots(subplot_kw=dict(xlabel='Time [s]', ylabel='%s field strength [V/m]'), num='rx' + str(rx), figsize=(20, 10), facecolor='w', edgecolor='w')
				line = ax.plot(time, outputdata,'r', lw=2, label=outputtext)
				ax.set_xlim([0, np.amax(time)])
		        #ax.set_ylim([-15, 20])
				ax.grid()
		        
				if 'H' in output:
					plt.setp(line, color='g')
					plt.setp(ax, ylabel='%s, field strength [A/m]' % outputtext)
				elif 'I' in output:
					plt.setp(line, color='b')
					plt.setp(ax, ylabel='%s, current [A]'% outputtext)

		# If multiple outputs required, create all nine subplots and populate only the specified ones
		else:
			fig, ax = plt.subplots(subplot_kw=dict(xlabel='Time [s]'), num='rx' + str(rx), figsize=(20, 10), facecolor='w', edgecolor='w')
			gs = gridspec.GridSpec(3, 3, hspace=0.3, wspace=0.3)
			for output in args.outputs:
		        
		        # Check for polarity of output and if requested output is in file
				if output[0] == 'm':
					polarity = -1
					outputtext = '-' + output[1:]
					output = output[1:]
				else:
					polarity = 1
					outputtext = output
		        
		        # Check if requested output is in file
				if output not in availableoutputs:
					raise CmdInputError('Output(s) requested to plot: {}, but available output(s) for receiver {} in the file: {}'.format(', '.join(args.outputs), rx, ', '.join(availableoutputs)))
		        
				
				outputdata = f['%s%s' %(path, output)][:] * polarity
		        
				if output == 'Ex':
					ax = plt.subplot(gs[0, 0])
					ax.plot(time, outputdata,'r', lw=2, label=outputtext)
					ax.set_ylabel('%s, field strength [V/m]' % outputtext)
		        #ax.set_ylim([-15, 20])
				elif output == 'Ey':
					ax = plt.subplot(gs[1, 0])
					ax.plot(time, outputdata,'r', lw=2, label=outputtext)
					ax.set_ylabel('%s, field strength [V/m]' % outputtext)
		        #ax.set_ylim([-15, 20])
				elif output == 'Ez':
					ax = plt.subplot(gs[2, 0])
					ax.plot(time, outputdata,'r', lw=2, label=outputtext)
					ax.set_ylabel('%s, field strength [V/m]' % outputtext)
		        #ax.set_ylim([-15, 20])
				elif output == 'Hx':
					ax = plt.subplot(gs[0, 1])
					ax.plot(time, outputdata,'g', lw=2, label=outputtext)
					ax.set_ylabel('%s, field strength [A/m]' % outputtext)
		        #ax.set_ylim([-0.03, 0.03])
				elif output == 'Hy':
					ax = plt.subplot(gs[1, 1])
					ax.plot(time, outputdata,'g', lw=2, label=outputtext)
					ax.set_ylabel('%s, field strength [A/m]' % outputtext)
		        #ax.set_ylim([-0.03, 0.03])
				elif output == 'Hz':
					ax = plt.subplot(gs[2, 1])
					ax.plot(time, outputdata,'g', lw=2, label=outputtext)
					ax.set_ylabel('%s, field strength [A/m]' % outputtext)
		        #ax.set_ylim([-0.03, 0.03])
				elif output == 'Ix':
					ax = plt.subplot(gs[0, 2])
					ax.plot(time, outputdata,'b', lw=2, label=outputtext)
					ax.set_ylabel('%s, current [A]' % outputtext)
				elif output == 'Iy':
					ax = plt.subplot(gs[1, 2])
					ax.plot(time, outputdata,'b', lw=2, label=outputtext)
					ax.set_ylabel('%s, current [A]' % outputtext)
				elif output == 'Iz':
					ax = plt.subplot(gs[2, 2])
					ax.plot(time, outputdata,'b', lw=2, label=outputtext)
					ax.set_ylabel('%s, current [A]' % outputtext)
			for ax in fig.axes:
				ax.set_xlim([0, np.amax(time)])
				ax.grid()

		# Save a PDF/PNG of the figure
		#fig.savefig(os.path.splitext(os.path.abspath(file))[0] + '_rx' + str(rx) + '.pdf', dpi=None, format='pdf', bbox_inches='tight', pad_inches=0.1)
		#fig.savefig(os.path.splitext(os.path.abspath(file))[0] + '_rx' + str(rx) + '.png', dpi=150, format='png', bbox_inches='tight', pad_inches=0.1)

	plt.show()

# 4-ex
/**
 * @(#)VarScan.java
 *
 * Copyright (c) 2009-2013 Daniel C. Koboldt and Washington University in St. Louis
 *
 * COPYRIGHT
 */

package net.sf.varscan;

//Import required packages //

import java.io.*;
import java.util.*;
import java.text.*;


/**
 * A set of tools for variant detection in next-generation sequence data.
 *
 * @version	2.4
 *
 * @author Daniel C. Koboldt https://sourceforge.net/u/dkoboldt
 *
 * <BR>
 * <pre>
 * COMMANDS
 * pileup2snp [pileup file] OPTIONS
 * 			Call SNPs from a pileup file that meet certain cutoffs
 * 			Input: 	Pileup file and parameters
 * 			Output: SNPs file with read counts and p-value
 *
 * pileup2indel [pileup file] OPTIONS
 * 			Call indels from a pileup file that meet certain cutoffs
 * 			Input: 	Pileup file and parameters
 * 			Output: Indels file with read counts and p-value
 *
 * pileup2cns [pileup file] OPTIONS
 * 			Call consensus genotypes (reference or variant) at sites with sufficient coverage
 * 			Input: 	Pileup file and parameters
 * 			Output: Consensus file with genotypes, read counts and p-values
 *
 * mpileup2cns [pileup file] OPTIONS
 * 			Call consensus genotypes (reference or variant) across one or more samples
 * 			Input: 	SAMtools mpileup file and parameters
 * 			Output: Consensus file with genotypes, read counts and p-values, or VCF file
 *
 * somatic [normal_pileup] [tumor_pileup] [output] OPTIONS
 * 			Determine somatic status of SNPs from normal/tumor pileup for positions
 * 			Input:	Normal pileup, tumor pileup, and positions file
 * 			Output: SNPs file with read counts and somatic status
 *
 * readcounts [pileup] --variants-file [positions] --output-file [output]
 * 			Obtain read counts for each allele of variants from a pileup file
 * 			Input:	Variants file and pileupfile
 * 			Output: Variants file with read counts for each allele
 *
 * filter [variant file] OPTIONS
 * 			Filter a set of SNPs/indels based on coverage, reads, p-value, etc.
 * 			Input:	SNPs file with read counts and p-value
 * 			Output: Filtered SNPs file with read counts and p-value
 *
 * somaticFilter [somatic-status file] OPTIONS
 * 			Filter VarScan Somatic/Germline/LOH calls for clusters and proximal indels
 * 			Input:	VarScan output for SNPs or Indels (varscan.output.snp)
 * 			Output: Variants passing all filters (varscan.output.snp.filter)
 *
 * fpFilter [variant-file] [readcount-file] OPTIONS
 * 			Apply the false-positive filter to VarScan variant calls
 * 			Input:	VarScan output for SNPs or Indels (varscan.output.snp) and bam-readcount output
 * 			Output: Variants passing all filters (varscan.output.snp.fpfilter)
 *
 * processSomatic [somatic-status file] OPTIONS
 * 			Process VarScan output by somatic status and confidence
 * 			Input:	VarScan output for SNPs or Indels (varscan.output.snp)
 * 			Output: Variants by somatic status (varscan.output.snp.Somatic)
 *
 * copyCaller [copynumber file] OPTIONS
 * 			Process VarScan copynumber output to adjust for GC and make preliminary calls
 * 			Input:	VarScan copynumber output (varscan.output.copynumber)
 * 			Output: Normalized copy number with preliminary calls (varscan.output.copynumber.called)
 *
 * compare [file1] [file2] [type] [output] OPTIONS
 * 			Compares chromosome-position entries in two tab-delimited files
 * 			Input:	File 1 and File 2
 * 			Output: Merged, intersected, or unique entries
 *
 * limit [variants] --regions-file [regions] --output-file [output]
 * 			Limit a tab-delimited file (SNPs, pileup, etc) to a set of positions or regions
 * 			Input:	tab-delimited input file with chromosome & position; positions-file or regions-file
 * 			Output: Entries in input-file matching regions or positions
 *
 * coverage [pileup-file] --regions-file [regions] --output-file [output]
 * 			**Experimental** Calculate Q>20 coverage depth/breadth for a set of target regions
 * 			Input:	Pileup file and tab-delimited regions-file
 * 			Output: Coverage report at various Q>20 depths (1x,10x,20x...)

 *
 * </pre>
 *
 *
 */
public class VarScan {

	static public class Copynumber {
		static public class SmartFileReader extends FileReader {

			public SmartFileReader(File file) throws FileNotFoundException {
				super(file);
			}

			public SmartFileReader(FileDescriptor fd) {
				super(fd);
			}

			public SmartFileReader(String fileName) throws FileNotFoundException {
				super(fileName);
			}

			public boolean ready() throws IOException {
				super.ready();
				// Implemented because sun.nio.cs.StreamDecoder doesn't implement ready() properly.
				return true;
			}

		}
		////////////////////////////////////////////////////////////////////////////////////////////////////
		// Constructor with two arguments (string[], boolean) expects mpileup input 					  //
		////////////////////////////////////////////////////////////////////////////////////////////////////
		public Copynumber(String[] args, boolean isMpileup)
		{
			String usage = "USAGE: java -jar VarScan.jar copynumber [normal-tumor.mpileup] [Opt: output] OPTIONS\n" +
					"\tnormal-tumor.mpileup - The SAMtools mpileup file for Normal and Tumor\n" +
					"\toutput - Output base name for files\n" +
					"\nOPTIONS:\n" +
					"\t--min-base-qual - Minimum base quality to count for coverage [20]\n" +
					"\t--min-map-qual - Minimum read mapping quality to count for coverage [20]\n" +
					"\t--min-coverage - Minimum coverage threshold for copynumber segments [20]\n" +
					"\t--min-segment-size - Minimum number of consecutive bases to report a segment [10]\n" +
					"\t--max-segment-size - Max size before a new segment is made [100]\n" +
					"\t--p-value - P-value threshold for significant copynumber change-point [0.01]\n" +
					"\t--data-ratio - The normal/tumor input data ratio for copynumber adjustment [1.0]\n";

			if(args.length < 2)
			{
				System.err.println(usage);
				return;
			}

			// Set parameter defaults //

			HashMap<String, String> params = VarScan.getParams(args);

			// Set up formatting for p-values //
			DecimalFormat pvalueFormat = new DecimalFormat("0.####E0");

			String outputName = "output";

			if(args.length >= 3 && !args[2].startsWith("-"))
			{
				outputName = args[2];
			}

			//	Set parameter defaults //

			int minCoverage = 	10;
			int minBaseQual = 	15;
			int minSegmentSize = 10;
			int maxSegmentSize = 100;
			double dataRatio = 1.00;
			double pValueThreshold = 0.01;
			long numBases = 0;

			// Try adjusting any provided parameters based on user inut //
			try
			{
				if(params.containsKey("min-coverage"))
				{
					minCoverage = Integer.parseInt(params.get("min-coverage"));
				}

				if(params.containsKey("min-base-qual"))
					minBaseQual = Integer.parseInt(params.get("min-base-qual"));

				if(params.containsKey("min-segment-size"))
					minSegmentSize = Integer.parseInt(params.get("min-segment-size"));

				if(params.containsKey("max-segment-size"))
					maxSegmentSize = Integer.parseInt(params.get("max-segment-size"));

				if(params.containsKey("p-value"))
					pValueThreshold = Double.parseDouble(params.get("p-value"));

				if(params.containsKey("data-ratio"))
					dataRatio = Double.parseDouble(params.get("data-ratio"));

				System.err.println("Min coverage:\t" + minCoverage);
				System.err.println("Min avg qual:\t" + minBaseQual);
				System.err.println("P-value thresh:\t" + pValueThreshold);

			}
			catch(Exception e)
			{
				System.err.println("Input Parameter Threw Exception: " + e.getLocalizedMessage());
				e.printStackTrace(System.err);
				System.exit(1);
			}

			// Print usage if -h or --help invoked //
			if(params.containsKey("help") || params.containsKey("h"))
			{
				System.err.println(usage);
				return;
			}

			// Check for correct input //

			if(args.length < 3)
			{
				System.err.println("Please provide an output file basename!");
				System.err.println(usage);
				System.exit(1);
			}


			// Parse piped input or user-provided pileup file //

			try
			{
				// Declare file-parsing variables //

				BufferedReader in = VarScan.getInfile(args);
				String line;

				// If no input, print usage //

				if(in == null)
				{
					System.out.println(usage);
					return;
				}

				// If input file not ready, give it a few seconds //
				int numNaps = 0;

				while(!in.ready())
				{
					try {
						Thread.sleep(5000);
						numNaps++;

						if(numNaps > 100)
						{
							System.err.println("Input file was not ready after 100 5-second cycles!");
							System.exit(10);
						}
					}
					catch(Exception e)
					{
						System.err.println("Exception while trying to get input" + e.getMessage());
						System.exit(1);
					}
				}

				// Proceed if input stream is ready //

				if(in != null && in.ready())
				{
					// Declare output file //
					PrintStream outCopySegments = null; // declare a print stream object for copynumber segments

					outCopySegments = new PrintStream( new FileOutputStream(outputName + ".copynumber") );
					outCopySegments.println("chrom\tchr_start\tchr_stop\tnum_positions\tnormal_depth\ttumor_depth\tlog2_ratio\tgc_content");


					System.err.println("Reading mpileup input...");
					int numParsingExceptions = 0;

					// Statistics counters //
					long sharedPositions = 0;
					long comparedPositions = 0;
					long rawCopySegments = 0;
					long goodCopySegments = 0;

					// Set some default parsing variables //
					String chromNormal = "";
					String chromTumor = "";
					String refBase = "";
					int posNormal = 0;
					int posTumor = 0;

					// Parameters for copy number calling //
					String copyChrom = "";
					int copyStart = 0;
					int copyStop = 0;
					int copyDepthNormal = 0;
					int copyDepthTumor = 0;
					long copySumNormal = 0;
					long copySumTumor = 0;
					long copyPositions = 0;
					long copyPositionsGC = 0;

					DecimalFormat oneDigit = new DecimalFormat("#0.0");
					DecimalFormat threeDigits = new DecimalFormat("#0.000");

					// Parse the infile line by line //

					while ((line = in.readLine()) != null)
					{
						numBases++;

						// Begin try-catch for line parsing //

						try
						{
							String[] lineContents = line.split("\t");

							// Verify expected pileup format //

//	    				if(lineContents.length > 5 && lineContents[0].length() > 0 && lineContents[1].length() > 0 && lineContents[2].length() > 0 && lineContents[3].length() > 0)
							if(lineContents.length < 8)
							{
								// This is an incomplete mpileup line, so skip it. If verbose, throw a warning //
								if(params.containsKey("verbose"))
								{
									System.err.println("Incomplete mpileup at line " + numBases + "; line being skipped.");
								}
							}
							else
							{
								sharedPositions++;

								// Parse common fields from line //
								String refName = lineContents[0];
								int position = Integer.parseInt(lineContents[1]);
								refBase = lineContents[2].toUpperCase();

								chromNormal = refName;
								chromTumor = refName;
								posNormal = position;
								posTumor = position;

								// Parse normal, which should be first sample //
								int normalOffset = 3;
								int pileupDepthNormal = 0;
								String normalQualities = "";
								if(lineContents.length >= (normalOffset + 1))
								{
									pileupDepthNormal = Integer.parseInt(lineContents[normalOffset]);
									//String normalBases = lineContents[normalOffset + 1];
									normalQualities = lineContents[normalOffset + 2];
								}

								// Parse tumor, which should be second sample //
								int tumorOffset = 6;
								int pileupDepthTumor = 0;
								String tumorQualities = "";
								if(lineContents.length >= (tumorOffset + 2 + 1))
								{
									pileupDepthTumor = Integer.parseInt(lineContents[tumorOffset]);
									//String tumorBases = lineContents[tumorOffset + 1];
									tumorQualities = lineContents[tumorOffset + 2];
								}


								// If either sample met the minimum coverage and both had at least one read //

//		    	        	if((pileupDepthNormal >= minCoverage || pileupDepthTumor >= minCoverage) && normalQualities.length() > 0)// && tumorQualities.length() > 0)

								// We want the normal sample to meet the minimum coverage because that's the comparator //
								if(pileupDepthNormal >= minCoverage && normalQualities.length() > 0)// && tumorQualities.length() > 0)
								{
									comparedPositions++;
									// Get the depth of bases above minimum quality //

									int normalDepth = VarScan.qualityDepth(normalQualities, minBaseQual);
									int tumorDepth = 0;
									if(tumorQualities.length() > 0)
										tumorDepth = VarScan.qualityDepth(tumorQualities, minBaseQual);

									// Determine if we have a copy changepoint //
									// If this base is not contiguous with the copyRegion
									// If the normal or tumor depth changes //

									int diffNormal = Math.abs(copyDepthNormal - normalDepth);
									int diffTumor = Math.abs(copyDepthTumor - tumorDepth);
									int posDiff = posTumor - copyStop;

									// DETERMINE IF WE CONTINUE THIS REGION OR PROCESS IT AND START A NEW ONE //

									boolean continueFlag = false;

									// If chromosomes differ or contiguity broken, process the region //

									if(posDiff > 2 || !(copyChrom.equals(chromTumor)))
									{
										continueFlag = false;
									}
									else
									{
										if(copyPositions >= maxSegmentSize)
										{
											continueFlag = false;
										}
										else if(diffNormal <= 2 && diffTumor <= 2)
										{
											continueFlag = true;
										}
										else
										{
											// Do a Fisher's exact test on the copy number changes. ##

											double changePvalue = VarScan.getSignificance(copyDepthNormal, copyDepthTumor, normalDepth, tumorDepth);

											// If depth change not significant, continue with region //
											if(changePvalue >= pValueThreshold)
											{
												continueFlag = true;
											}
											else
											{
												continueFlag = false;
											}

										}
									}


									// If continuing, extend this region and don't process yet //

									if(continueFlag)
									{
										copySumNormal += normalDepth;
										copySumTumor += tumorDepth;
										copyPositions++;
										if(refBase.equals("G") || refBase.equals("C") || refBase.equals("g") || refBase.equals("c"))
											copyPositionsGC++;
										copyStop = posTumor;
									}

									// Otherwise, process this region (if it qualifies) and start a new one //

									else
									{
										if(copyPositions >= minSegmentSize)
										{
											rawCopySegments++;
											String regionResults = processCopyRegion(copyChrom, copyStart, copyStop, copyPositions, copyPositionsGC, copySumNormal, copySumTumor, minCoverage, dataRatio);

											if(regionResults.length() > 0)
											{
												outCopySegments.println(regionResults);
												goodCopySegments++;
											}
										}

										// Start a new copyNumber region //
										copyChrom = chromTumor;
										copyStart = posTumor;
										copyStop = posTumor;
										copyDepthNormal = normalDepth;
										copyDepthTumor = tumorDepth;
										copySumNormal = normalDepth;
										copySumTumor = tumorDepth;
										copyPositions = 1;
										if(refBase.equals("G") || refBase.equals("C") || refBase.equals("g") || refBase.equals("c"))
											copyPositionsGC = 1;
										else
											copyPositionsGC = 0;
									}


								}
								else
								{
									// If minimum coverage was not met, print region //
									// If we had a copyNumber region that met minimum coverage, report it //
									if(copyPositions >= minSegmentSize)
									{
										rawCopySegments++;
										String regionResults = processCopyRegion(copyChrom, copyStart, copyStop, copyPositions, copyPositionsGC, copySumNormal, copySumTumor, minCoverage, dataRatio);

										if(regionResults.length() > 0)
										{
											outCopySegments.println(regionResults);
											goodCopySegments++;
										}
									}

									// Reset the copyNumber region //
									copyChrom = "";
									copyStart = 0;
									copyStop = 0;
									copyDepthNormal = 0;
									copyDepthTumor = 0;
									copySumNormal = 0;
									copySumTumor = 0;
									copyPositions = 0;
									copyPositionsGC = 0;
								}

							}

						}
						catch(Exception e)
						{
							System.err.println("Parsing Exception on line:\n" + line + "\n" + e.getLocalizedMessage());
							numParsingExceptions++;
							if(numParsingExceptions >= 5)
							{
								System.err.println("Too many parsing exceptions encountered; exiting");
								return;
							}
							return;
						}


					}

					// Last region: If minimum coverage was not met, print region //
					// If we had a copyNumber region that met minimum coverage, report it //
					if(copyPositions > minSegmentSize)
					{
						rawCopySegments++;
						String regionResults = processCopyRegion(copyChrom, copyStart, copyStop, copyPositions, copyPositionsGC, copySumNormal, copySumTumor, minCoverage, dataRatio);

						if(regionResults.length() > 0)
						{
							outCopySegments.println(regionResults);
							goodCopySegments++;
						}
					}

					in.close();

					System.err.println(sharedPositions + " positions in mpileup"); //stats.get("sharedPositions")
					System.err.println(comparedPositions + " had sufficient coverage for comparison"); //stats.get("comparedPositions")
					System.err.println(rawCopySegments + " raw copynumber segments with size > " + minSegmentSize);
					System.err.println(goodCopySegments + " good copynumber segments with depth > " + minCoverage);

				}
				else
				{
					System.err.println("Input file never ready for parsing (maybe due to file I/O)...");
					System.exit(10);
				}
			}
			catch (IOException e)
			{
				System.err.println("File Parsing Exception: " + e.getLocalizedMessage());
				e.printStackTrace(System.err);
				System.exit(11);
			}



		}



		////////////////////////////////////////////////////////////////////////////////////////////////////
		// Constructor with one argument (string[]) expects independent normal and tumor pileups as input //
		////////////////////////////////////////////////////////////////////////////////////////////////////

		public Copynumber(String[] args)
		{
			String usage = "USAGE: VarScan copynumber [normal_pileup] [tumor_pileup] [Opt: output] OPTIONS\n" +
					"\tnormal_pileup - The SAMtools pileup file for Normal\n" +
					"\ttumor_pileup - The SAMtools pileup file for Tumor\n" +
					"\toutput - Output base name for files\n" +
					"***If you have a single mpileup, see VarScan copynumber -mpileup 1 -h ***\n" +
					"\nOPTIONS:\n" +
					"\t--min-base-qual - Minimum base quality to count for coverage [20]\n" +
					"\t--min-map-qual - Minimum read mapping quality to count for coverage [20]\n" +
					"\t--min-coverage - Minimum coverage threshold for copynumber segments [20]\n" +
					"\t--min-segment-size - Minimum number of consecutive bases to report a segment [10]\n" +
					"\t--max-segment-size - Max size before a new segment is made [100]\n" +
					"\t--p-value - P-value threshold for significant copynumber change-point [0.01]\n" +
					"\t--data-ratio - The normal/tumor input data ratio for copynumber adjustment [1.0]\n";

			if(args.length < 3)
			{
				System.err.println(usage);
				return;
			}

			// Get the required arguments //
			String normalPileupFile = args[1];
			String tumorPileupFile = args[2];

			String outputName = "output";

			if(args.length >= 4 && !args[3].startsWith("-"))
			{
				outputName = args[3];
			}

			System.err.println("Normal Pileup: " + normalPileupFile);
			System.err.println("Tumor Pileup: " + tumorPileupFile);
			System.err.println("NOTICE: While dual input files are still supported, using a single mpileup file (normal-tumor) with the --mpileup 1 setting is strongly recommended.");

			//	Set parameter defaults //

			int minCoverage = 	10;
			int minBaseQual = 	15;
			int minSegmentSize = 10;
			int maxSegmentSize = 100;
			double dataRatio = 1.00;
			double pValueThreshold = 0.01;

			// Parse command-line parameters //
			HashMap<String, String> params = VarScan.getParams(args);

			// Try adjusting any provided parameters based on user inut //
			try
			{
				if(params.containsKey("min-coverage"))
				{
					minCoverage = Integer.parseInt(params.get("min-coverage"));
				}

				if(params.containsKey("min-base-qual"))
					minBaseQual = Integer.parseInt(params.get("min-base-qual"));

				if(params.containsKey("min-segment-size"))
					minSegmentSize = Integer.parseInt(params.get("min-segment-size"));

				if(params.containsKey("max-segment-size"))
					maxSegmentSize = Integer.parseInt(params.get("max-segment-size"));

				if(params.containsKey("p-value"))
					pValueThreshold = Double.parseDouble(params.get("p-value"));

				if(params.containsKey("data-ratio"))
					dataRatio = Double.parseDouble(params.get("data-ratio"));

				System.err.println("Min coverage:\t" + minCoverage);
				System.err.println("Min avg qual:\t" + minBaseQual);
				System.err.println("P-value thresh:\t" + pValueThreshold);

			}
			catch(Exception e)
			{
				System.err.println("Input Parameter Threw Exception: " + e.getLocalizedMessage());
				e.printStackTrace(System.err);
				System.exit(1);
			}

			// Print usage if -h or --help invoked //
			if(params.containsKey("help") || params.containsKey("h"))
			{
				System.err.println(usage);
				return;
			}

			// Check for correct input //

			if(args.length < 3)
			{
				System.err.println("Please provide an output file basename!");
				System.err.println(usage);
				System.exit(1);
			}


			// Statistics counters //
			long tumorPositions = 0;
			long sharedPositions = 0;
			long comparedPositions = 0;
			long rawCopySegments = 0;
			long goodCopySegments = 0;

			try
			{
				// Declare output file //
				PrintStream outCopySegments = null; // declare a print stream object for copynumber segments

				outCopySegments = new PrintStream( new FileOutputStream(outputName + ".copynumber") );
				outCopySegments.println("chrom\tchr_start\tchr_stop\tnum_positions\tnormal_depth\ttumor_depth\tlog2_ratio\tgc_content");

				// Prepare file readers for normal and tumor pileups //

				BufferedReader normal = new BufferedReader(new SmartFileReader(normalPileupFile));
				BufferedReader tumor = new BufferedReader(new SmartFileReader(tumorPileupFile));

				if(!(normal.ready() && tumor.ready()))
				{
					// Delay a few seconds to let SAMtools pileup start outputting //
					try {
						Thread.sleep(5000);

						if(!(normal.ready() && tumor.ready()))
							Thread.sleep(5000);

						if(!(normal.ready() && tumor.ready()))
							Thread.sleep(5000);

						if(!(normal.ready() && tumor.ready()))
							Thread.sleep(5000);
					}
					catch(Exception e)
					{

					}

				}

				// Exit if files not ready after waiting //

				if(!(normal.ready() && tumor.ready()))
				{
					System.err.println("ERROR: Invalid input file(s)");
					System.exit(10);
				}

				String lineNormal;
				String lineTumor;
				String chromNormal = "";
				String chromTumor = "";
				String prevChromNormal = "";
				String prevChromTumor = "";
				String refBase = "";
				int posNormal = 0;
				int posTumor = 0;

				// Parameters for copy number calling //
				String copyChrom = "";
				int copyStart = 0;
				int copyStop = 0;
				int copyDepthNormal = 0;
				int copyDepthTumor = 0;
				long copySumNormal = 0;
				long copySumTumor = 0;
				long copyPositions = 0;
				long copyPositionsGC = 0;

				DecimalFormat oneDigit = new DecimalFormat("#0.0");
				DecimalFormat threeDigits = new DecimalFormat("#0.000");

				try {
					// Get first line of Normal //

					if((lineNormal = normal.readLine()) != null)
					{
						String[] normalContents = lineNormal.split("\t");

						if(normalContents.length > 1)
						{
							chromNormal = normalContents[0];
							posNormal = Integer.parseInt(normalContents[1]);
						}
					}

					// Loop through lines in tumor //

					while ((lineTumor = tumor.readLine()) != null)
					{
						tumorPositions++;
						String[] tumorContents = lineTumor.split("\t");

						if(tumorContents.length > 1)
						{
							chromTumor = tumorContents[0];
							posTumor = Integer.parseInt(tumorContents[1]);
						}

						// Parse normal lines until we get the same chromosome //
						boolean flagEOF = false;
						boolean normalWasReset = false;

						//	Advance in normal file if tumor is changed but normal is not, or if tumor is higher //
						while(!chromNormal.equals(chromTumor) && !chromTumor.equals(prevChromTumor) && !flagEOF && (chromNormal.equals(prevChromTumor) || inSortOrder(chromNormal, chromTumor)))
						{
							//System.err.println("Normal (" + chromNormal + ") catching up to " + chromTumor);
							// Get next line from normal pileup //
							if((lineNormal = normal.readLine()) != null)
							{
								String[] normalContents = lineNormal.split("\t");

								if(normalContents.length > 1)
								{
									chromNormal = normalContents[0];
									posNormal = Integer.parseInt(normalContents[1]);
								}
							}
							else
							{
								flagEOF = true;
							}


						}

						// If chromosomes match and are non-blank, attempt to get matching positions //
						if(chromNormal.equals(chromTumor) && !chromNormal.equals(""))
						{
							normalWasReset = false;
							// Seek to matching Normal Position //

							while(chromNormal.equals(chromTumor) && posNormal < posTumor && ((lineNormal = normal.readLine()) != null))
							{
								String[] normalContents = lineNormal.split("\t");
								if(normalContents.length > 1)
								{
									chromNormal = normalContents[0];
									posNormal = Integer.parseInt(normalContents[1]);

									// If still less than tumor position, look for homozygous del //
									if(posNormal < posTumor)
									{
										int pileupDepthNormal = 0;
										String normalQualities = "";

										try
										{
											// Pileup Files have 6-7 columns //
											if(normalContents.length <= 7)
											{
												pileupDepthNormal = Integer.parseInt(normalContents[3]);
												normalQualities = normalContents[5];
											}
											// Pileup lines in CNS files have 10-11 columns
											else if (normalContents.length >= 10 && normalContents.length <= 11)
											{
												pileupDepthNormal = Integer.parseInt(normalContents[7]);
												normalQualities = normalContents[9];
											}
										}
										catch(Exception e)
										{

										}

									}
									else
									{

									}
								}
							}

							// Seek to matching Tumor Position //

							while(chromNormal.equals(chromTumor) && posTumor < posNormal && ((lineTumor = tumor.readLine()) != null))
							{
								tumorContents = lineTumor.split("\t");
								if(tumorContents.length > 1)
								{
									chromTumor = tumorContents[0];
									posTumor = Integer.parseInt(tumorContents[1]);
								}
							}

							// Proceed if normal and tumor positions match //

							if(chromNormal.equals(chromTumor) && chromNormal.equals(chromTumor) && posNormal == posTumor)
							{
								//stats.put("sharedPositions", (stats.get("sharedPositions") + 1));
								sharedPositions++;
								refBase = tumorContents[2];

//			    			 Parse out base qualities //
								String[] normalContents = lineNormal.split("\t");
								int pileupDepthNormal = 0;
								int pileupDepthTumor = 0;
								String normalQualities = "";
								String tumorQualities = "";

								// Pileup Files have 6-7 columns //
								if(normalContents.length >= 6 && normalContents.length <= 7)
								{
									pileupDepthNormal = Integer.parseInt(normalContents[3]);
									normalQualities = normalContents[5];
								}
								// Pileup lines in CNS files have 10-11 columns
								else if (normalContents.length >= 10 && normalContents.length <= 11)
								{
									pileupDepthNormal = Integer.parseInt(normalContents[7]);
									normalQualities = normalContents[9];
								}

								// Pileup Files have 6-7 columns //
								if(tumorContents.length >= 6 && tumorContents.length <= 7)
								{
									tumorQualities = tumorContents[5];
									pileupDepthTumor = Integer.parseInt(tumorContents[3]);
								}
								// Pileup lines in CNS files have 10-11 columns
								else if (tumorContents.length >= 10 && tumorContents.length <= 11)
								{
									tumorQualities = tumorContents[9];
									pileupDepthTumor = Integer.parseInt(tumorContents[7]);
								}

								// If either sample met the minimum coverage and both had at least one read //

//	    					if((pileupDepthNormal >= minCoverage || pileupDepthTumor >= minCoverage) && normalQualities.length() > 0 && tumorQualities.length() > 0)

								// We want the normal sample to meet the minimum coverage because that's the comparator //
								if(pileupDepthNormal >= minCoverage && normalQualities.length() > 0) // && tumorQualities.length() > 0)
								{
									comparedPositions++;
//	    						 Get the depth of bases above minimum quality //

									int normalDepth = VarScan.qualityDepth(normalQualities, minBaseQual);
									int tumorDepth = VarScan.qualityDepth(tumorQualities, minBaseQual);

									// Determine if we have a copy changepoint //
									// If this base is not contiguous with the copyRegion
									// If the normal or tumor depth changes //

									int diffNormal = Math.abs(copyDepthNormal - normalDepth);
									int diffTumor = Math.abs(copyDepthTumor - tumorDepth);
									int posDiff = posTumor - copyStop;

									// DETERMINE IF WE CONTINUE THIS REGION OR PROCESS IT AND START A NEW ONE //

									boolean continueFlag = false;

									// If chromosomes differ or contiguity broken, process the region //

									if(posDiff > 2 || !(copyChrom.equals(chromTumor)))
									{
										continueFlag = false;
									}
									else
									{
										if(copyPositions >= maxSegmentSize)
										{
											continueFlag = false;
										}
										else if(diffNormal <= 2 && diffTumor <= 2)
										{
											continueFlag = true;
										}
										else
										{
											// Do a Fisher's exact test on the copy number changes. ##

											double changePvalue = VarScan.getSignificance(copyDepthNormal, copyDepthTumor, normalDepth, tumorDepth);

											// If depth change not significant, continue with region //
											if(changePvalue >= pValueThreshold)
											{
												continueFlag = true;
											}
											else
											{
												continueFlag = false;
											}

										}
									}


									// If continuing, extend this region and don't process yet //

									if(continueFlag)
									{
										copySumNormal += normalDepth;
										copySumTumor += tumorDepth;
										copyPositions++;
										if(refBase.equals("G") || refBase.equals("C") || refBase.equals("g") || refBase.equals("c"))
											copyPositionsGC++;
										copyStop = posTumor;
									}

									// Otherwise, process this region (if it qualifies) and start a new one //

									else
									{
										if(copyPositions >= minSegmentSize)
										{
											rawCopySegments++;
											String regionResults = processCopyRegion(copyChrom, copyStart, copyStop, copyPositions, copyPositionsGC, copySumNormal, copySumTumor, minCoverage, dataRatio);

											if(regionResults.length() > 0)
											{
												outCopySegments.println(regionResults);
												goodCopySegments++;
											}
										}

										// Start a new copyNumber region //
										copyChrom = chromTumor;
										copyStart = posTumor;
										copyStop = posTumor;
										copyDepthNormal = normalDepth;
										copyDepthTumor = tumorDepth;
										copySumNormal = normalDepth;
										copySumTumor = tumorDepth;
										copyPositions = 1;
										if(refBase.equals("G") || refBase.equals("C") || refBase.equals("g") || refBase.equals("c"))
											copyPositionsGC = 1;
										else
											copyPositionsGC = 0;
									}


								}
								else
								{
									// If minimum coverage was not met, print region //
									// If we had a copyNumber region that met minimum coverage, report it //
									if(copyPositions >= minSegmentSize)
									{
										rawCopySegments++;
										String regionResults = processCopyRegion(copyChrom, copyStart, copyStop, copyPositions, copyPositionsGC, copySumNormal, copySumTumor, minCoverage, dataRatio);

										if(regionResults.length() > 0)
										{
											outCopySegments.println(regionResults);
											goodCopySegments++;
										}
									}

									// Reset the copyNumber region //
									copyChrom = "";
									copyStart = 0;
									copyStop = 0;
									copyDepthNormal = 0;
									copyDepthTumor = 0;
									copySumNormal = 0;
									copySumTumor = 0;
									copyPositions = 0;
									copyPositionsGC = 0;
								}

								// Record this chromosome //

								prevChromNormal = chromNormal;
								prevChromTumor = chromTumor;
							}
							else
							{
								//System.err.println("Failed to match positions " + chromNormal + " " + posNormal + " to Tumor " + chromTumor + " " + posTumor);
							}
						}
						// If they're in sort order, do nothing so that tumor can catch up //
						else if(inSortOrder(chromNormal, chromTumor))
						{
							System.err.println("Not resetting normal file because " + chromNormal + " < " + chromTumor);
						}
						// If we reached the end of the normal file but never saw this chromosome, //
						// fast-forward until tumor chromosome changes and reset normal file //
						else if(flagEOF)
						{
							flagEOF = false;

							while(prevChromTumor.equals(chromTumor) && !flagEOF)
							{
								if((lineTumor = tumor.readLine()) != null)
								{
									tumorContents = lineTumor.split("\t");

									if(tumorContents.length > 1)
									{
										chromTumor = tumorContents[0];
										posTumor = Integer.parseInt(tumorContents[1]);
									}
								}
								else
								{
									flagEOF = true;
								}
							}

							// Reset the normal file if we've already passed this chromosome in normal //

							if(!flagEOF && !normalWasReset)
							{
								if(inSortOrder(chromNormal, chromTumor))
								{
									System.err.println("Not resetting normal file because " + chromNormal + " < " + chromTumor);
								}
								else
								{
									System.err.println("Resetting normal file because " + chromNormal + " > " + chromTumor);
									normalWasReset = true;
									normal.close();
									normal = new BufferedReader(new SmartFileReader(normalPileupFile));
								}

							}
						}

					}

				}
				catch (Exception e)
				{
					System.err.println("Exception encountered while parsing normal/tumor files: " + e.getMessage());
					System.err.println("Note: It is HIGHLY recommended that you use a two-sample mpileup input rather than separate pileup files for normal/tumor.");
					return;
				}





				normal.close();
				tumor.close();

				// If we had a copyNumber region that met minimum coverage, report it //
				if(copyPositions > minSegmentSize)
				{
					rawCopySegments++;
					String regionResults = processCopyRegion(copyChrom, copyStart, copyStop, copyPositions, copyPositionsGC, copySumNormal, copySumTumor, minCoverage, dataRatio);

					if(regionResults.length() > 0)
					{
						outCopySegments.println(regionResults);
						goodCopySegments++;
					}
				}


				outCopySegments.close();

				System.err.println(tumorPositions + " positions in tumor");
				System.err.println(sharedPositions + " positions shared in normal"); //stats.get("sharedPositions")
				System.err.println(comparedPositions + " had sufficient coverage for comparison"); //stats.get("comparedPositions")

				System.err.println(rawCopySegments + " raw copynumber segments with size > " + minSegmentSize);
				System.err.println(goodCopySegments + " good copynumber segments with depth > " + minCoverage);
			}
			catch (IOException e)
			{
				System.err.println("File Parsing Exception: " + e.getLocalizedMessage());
				e.printStackTrace(System.err);
				System.exit(11);
			}
		}


		/**
		 * Calculates relative tumor copynumber for a contiguous segment
		 *
//		 * @param	args	Command-line arguments
		 * @return			HashMap of parameter names and their values
		 */
		static String processCopyRegion(String copyChrom, int copyStart, int copyStop, long copyPositions, long copyPositionsGC, long copySumNormal, long copySumTumor, int minCoverage, double dataRatio)
		{
			DecimalFormat oneDigit = new DecimalFormat("#0.0");
			DecimalFormat threeDigits = new DecimalFormat("#0.000");

			try
			{
//			 Calculate average depth //
				float avgNormal = (float) copySumNormal / (float) copyPositions;
				float avgTumor = (float) copySumTumor / (float) copyPositions;
				// Adjust tumor depth for ratio
				float adjustedTumorDepth = (float) dataRatio * (float) avgTumor;

				float gcContent = (float) copyPositionsGC / (float) copyPositions * 100;

				if(avgNormal >= minCoverage || avgTumor >= minCoverage)
				{
					// Determine ratio and diff //
					if(avgNormal >= 0.01 && avgTumor >= 0.01)
					{
						float tumorNormalRatio = adjustedTumorDepth / avgNormal;
						double log2ratio = Math.log(tumorNormalRatio) / Math.log(2);

						return(copyChrom + "\t" + copyStart + "\t" + copyStop + "\t" + copyPositions + "\t" + oneDigit.format(avgNormal) + "\t" + oneDigit.format(avgTumor) + "\t" + threeDigits.format(log2ratio) + "\t" + oneDigit.format(gcContent));
					}
					else if (avgTumor >= 0.01)
					{
						// If only tumor has coverage, handle it //
						double log2ratio = 2.00;
						return(copyChrom + "\t" + copyStart + "\t" + copyStop + "\t" + copyPositions + "\t" + oneDigit.format(avgNormal) + "\t" + oneDigit.format(avgTumor) + "\t" + threeDigits.format(log2ratio) + "\t" + oneDigit.format(gcContent));
					}
					else
					{
						// If only normal has coverage, mark as homozygyous deletion //
						double log2ratio = -2.00;
						return(copyChrom + "\t" + copyStart + "\t" + copyStop + "\t" + copyPositions + "\t" + oneDigit.format(avgNormal) + "\t" + oneDigit.format(avgTumor) + "\t" + threeDigits.format(log2ratio) + "\t" + oneDigit.format(gcContent));
					}

				}
				else
				{
//				System.err.println("Warning: Not reporting region " + copyChrom + " " + copyStart + " " + copyStop + " " + copyPositions + " " + avgNormal + " " + avgTumor);
				}
			}
			catch(Exception e)
			{
				System.err.println("Warning: Error while processing copynumber segment:" + e.getMessage());
			}


			return("");
		}


		/**
		 * Determine if tumor chromosome is before normal chromosome in sort order
		 *
//		 * @param	args	Command-line arguments
		 * @return			HashMap of parameter names and their values
		 */
		static boolean inSortOrder(String chrom1, String chrom2)
		{
			String[] testArray = {chrom1, chrom2};
			Arrays.sort(testArray);

			if(testArray[0].equals(chrom1))
				return true;

			return false;
		}



		/**
		 * Determines the sort order for chromosomes
		 *
//		 * @param	args	Command-line arguments
		 * @return			HashMap of parameter names and their values
		 */
		static Boolean chromSorted(String chrom1, String chrom2)
		{
			Boolean answer = false;

			chrom1.replace("X", "23");
			chrom1.replace("Y", "24");
			chrom1.replace("M", "25");

			chrom2.replace("X", "23");
			chrom2.replace("Y", "24");
			chrom2.replace("M", "25");

			String[] unsorted = {chrom1, chrom2};
			String[] sorted = {chrom1, chrom2};
			Arrays.sort(sorted);
			System.err.println("Sorted order is " + sorted[0] + " " + sorted[1]);
			try{
				if(sorted[0].equals(unsorted[0]))
				{
					answer = true;
				}
			}
			catch(Exception e)
			{

			}

			return(answer);
		}
	}

	public static class FishersExact {
		private static final boolean DEBUG = false;
		private double[] f;
		int maxSize;


		/**
		 * constructor for FisherExact table
		 *
		 * @param maxSize is the maximum sum that will be encountered by the table (a+b+c+d)
		 */
		public FishersExact(int maxSize) {
			this.maxSize = maxSize;
			double cf = 1.0;
			f = new double[maxSize + 1];
			f[0] = 0.0;
			for (int i = 1; i <= this.maxSize; i++) {
				f[i] = f[i - 1] + Math.log(i);
			}
		}

		/**
		 * calculates the P-value for this specific state
		 *
		 * @param a     a, b, c, d are the four cells in a 2x2 matrix
		 * @param b
		 * @param c
		 * @param d
		 * @return the P-value
		 */
		public final double getP(int a, int b, int c, int d) {
			try
			{
				int n = a + b + c + d;
				if (n > maxSize) {
					return Double.NaN;
				}
				double p;
				p = (f[a + b] + f[c + d] + f[a + c] + f[b + d]) - (f[a] + f[b] + f[c] + f[d] + f[n]);
				return Math.exp(p);
			}
			catch(Exception e)
			{
				return Double.NaN;
			}

		}

		/**
		 * Calculates the one-tail P-value for the Fisher Exact test.  Determines whether to calculate the right- or left-
		 * tail, thereby always returning the smallest p-value.
		 *
		 * @param a     a, b, c, d are the four cells in a 2x2 matrix
		 * @param b
		 * @param c
		 * @param d
		 * @return one-tailed P-value (right or left, whichever is smallest)
		 */
		public final double getCumlativeP(int a, int b, int c, int d) {
			int min, i;
			int n = a + b + c + d;
			if (n > maxSize) {
				return Double.NaN;
			}
			double p = 0;

			p += getP(a, b, c, d);
			if (DEBUG) {System.out.println("p = " + p);}
			if ((a * d) >= (b * c)) {
				if (DEBUG) {System.out.println("doing R-tail: a=" + a + " b=" + b + " c=" + c + " d=" + d);}
				min = (c < b) ? c : b;
				for (i = 0; i < min; i++) {
					if (DEBUG) {System.out.print("doing round " + i);}
					p += getP(++a, --b, --c, ++d);
					if (DEBUG) {System.out.println("\ta=" + a + " b=" + b + " c=" + c + " d=" + d);}
				}
				System.out.println("");
			}
			if ((a * d) < (b * c)) {
				if (DEBUG) {System.out.println("doing L-tail: a=" + a + " b=" + b + " c=" + c + " d=" + d);}
				min = (a < d) ? a : d;
				for (i = 0; i < min; i++) {
					if (DEBUG) {System.out.print("doing round " + i);}
					double pTemp = getP(--a, ++b, ++c, --d);
					if (DEBUG) {System.out.print("\tpTemp = " + pTemp);}
					p += pTemp;
					if (DEBUG) {System.out.println("\ta=" + a + " b=" + b + " c=" + c + " d=" + d);}
				}
			}
			return p;
		}

		/**
		 * Calculates the right-tail P-value for the Fisher Exact test.
		 *
		 * @param a     a, b, c, d are the four cells in a 2x2 matrix
		 * @param b
		 * @param c
		 * @param d
		 * @return one-tailed P-value (right-tail)
		 */
		public final double getRightTailedP(int a, int b, int c, int d) {
			int min, i;
			int n = a + b + c + d;
			if (n > maxSize) {
				return Double.NaN;
			}
			double p = 0;

			p += getP(a, b, c, d);
			if (DEBUG) {System.out.println("p = " + p);}
			if (DEBUG) {System.out.println("doing R-tail: a=" + a + " b=" + b + " c=" + c + " d=" + d);}
			min = (c < b) ? c : b;
			for (i = 0; i < min; i++) {
				p += getP(++a, --b, --c, ++d);

			}
			return p;
		}

		/**
		 * Calculates the left-tail P-value for the Fisher Exact test.
		 *
		 * @param a     a, b, c, d are the four cells in a 2x2 matrix
		 * @param b
		 * @param c
		 * @param d
		 * @return one-tailed P-value (left-tail)
		 */
		public final double getLeftTailedP(int a, int b, int c, int d) {
			int min, i;
			int n = a + b + c + d;
			if (n > maxSize) {
				return Double.NaN;
			}
			double p = 0;

			p += getP(a, b, c, d);
			if (DEBUG) {System.out.println("p = " + p);}
			if (DEBUG) {System.out.println("doing L-tail: a=" + a + " b=" + b + " c=" + c + " d=" + d);}
			min = (a < d) ? a : d;
			for (i = 0; i < min; i++) {
				if (DEBUG) {System.out.print("doing round " + i);}
				double pTemp = getP(--a, ++b, ++c, --d);
				if (DEBUG) {System.out.print("\tpTemp = " + pTemp);}
				p += pTemp;
				if (DEBUG) {System.out.println("\ta=" + a + " b=" + b + " c=" + c + " d=" + d);}
			}


			return p;
		}


		/**
		 *   Calculates the two-tailed P-value for the Fisher Exact test.
		 *
		 *   In order for a table under consideration to have its p-value included
		 *   in the final result, it must have a p-value less than the original table's P-value, i.e.
		 *   Fisher's exact test computes the probability, given the observed marginal
		 *   frequencies, of obtaining exactly the frequencies observed and any configuration more extreme.
		 *   By "more extreme," we mean any configuration (given observed marginals) with a smaller probability of
		 *   occurrence in the same direction (one-tailed) or in both directions (two-tailed).
		 *
		 * @param a     a, b, c, d are the four cells in a 2x2 matrix
		 * @param b
		 * @param c
		 * @param d
		 * @return two-tailed P-value
		 */
		public final double getTwoTailedP(int a, int b, int c, int d) {
			int min, i;
			int n = a + b + c + d;
			if (n > maxSize) {
				return Double.NaN;
			}
			double p = 0;

			double baseP = getP(a, b, c, d);
//         in order for a table under consideration to have its p-value included
//         in the final result, it must have a p-value less than the baseP, i.e.
//         Fisher's exact test computes the probability, given the observed marginal
//         frequencies, of obtaining exactly the frequencies observed and any configuration more extreme.
//         By "more extreme," we mean any configuration (given observed marginals) with a smaller probability of
//         occurrence in the same direction (one-tailed) or in both directions (two-tailed).

			if (DEBUG) {System.out.println("baseP = " + baseP);}
			int initialA = a, initialB = b, initialC = c, initialD = d;
			p += baseP;
			if (DEBUG) {System.out.println("p = " + p);}
			if (DEBUG) {System.out.println("Starting with R-tail: a=" + a + " b=" + b + " c=" + c + " d=" + d);}
			min = (c < b) ? c : b;
			for (i = 0; i < min; i++) {
				if (DEBUG) {System.out.print("doing round " + i);}
				double tempP = getP(++a, --b, --c, ++d);
				if (tempP <= baseP) {
					if (DEBUG) {System.out.print("\ttempP (" + tempP + ") is less than baseP (" + baseP + ")");}
					p += tempP;
				}
				if (DEBUG) {System.out.println(" a=" + a + " b=" + b + " c=" + c + " d=" + d);}
			}

			// reset the values to their original so we can repeat this process for the other side
			a = initialA;
			b = initialB;
			c = initialC;
			d = initialD;

			if (DEBUG) {System.out.println("Now doing L-tail: a=" + a + " b=" + b + " c=" + c + " d=" + d);}
			min = (a < d) ? a : d;
			if (DEBUG) {System.out.println("min = " + min);}
			for (i = 0; i < min; i++) {
				if (DEBUG) {System.out.print("doing round " + i);}
				double pTemp = getP(--a, ++b, ++c, --d);
				if (DEBUG) {System.out.println("  pTemp = " + pTemp);}
				if (pTemp <= baseP) {
					if (DEBUG) {System.out.print("\ttempP (" + pTemp + ") is less than baseP (" + baseP + ")");}
					p += pTemp;
				}
				if (DEBUG) {System.out.println(" a=" + a + " b=" + b + " c=" + c + " d=" + d);}
			}
			return p;
		}

//		public static void main(String[] args) {
//
//			int[][] argInts = new int[15][4];
//			argInts[0] = new int[]{2, 3, 6, 4};
//			argInts[1] = new int[]{2, 1, 3, 0};
//			argInts[2] = new int[]{3, 0, 2, 1};
//			argInts[3] = new int[]{1, 2, 0, 3};
//			argInts[4] = new int[]{3, 1, 1, 3};
//			argInts[5] = new int[]{1, 3, 3, 1};
//			argInts[6] = new int[]{0, 1, 1, 0};
//			argInts[7] = new int[]{1, 0, 0, 1};
//			argInts[8] = new int[]{11, 0, 0, 6};
//			argInts[9] = new int[]{10, 1, 1, 5};
//			argInts[10] = new int[]{5, 6, 6, 0};
//			argInts[11] = new int[]{9, 2, 2, 4};
//			argInts[12] = new int[]{6, 5, 5, 1};
//			argInts[13] = new int[]{8, 3, 3, 3};
//			argInts[14] = new int[]{7, 4, 4, 2};
//
//			net.sf.varscan.FishersExact fe = new net.sf.varscan.FishersExact(100);
//
//			for (int i = 0; i < argInts.length; i++) {
//				System.out.println("\na=" + argInts[i][0] + " b=" + argInts[i][1] + " c=" + argInts[i][2] + " d=" + argInts[i][3]);
//				System.out.print("*****Original algorithm: ");
//				double cumulativeP = fe.getCumlativeP(argInts[i][0], argInts[i][1], argInts[i][2], argInts[i][3]);
//				System.out.println("\tcumulativeP = " + cumulativeP);
//
//				System.out.print("*****Left Tailed: ");
//				double leftTailedP = fe.getLeftTailedP(argInts[i][0], argInts[i][1], argInts[i][2], argInts[i][3]);
//				System.out.println("\tleftTailedP = " + leftTailedP);
//
//				System.out.print("*****Right Tailed: ");
//				double rightTailedP = fe.getRightTailedP(argInts[i][0], argInts[i][1], argInts[i][2], argInts[i][3]);
//				System.out.println("\trightTailedP = " + rightTailedP);
//
//				System.out.print("*****Two Tailed: ");
//				double twoTailedP = fe.getTwoTailedP(argInts[i][0], argInts[i][1], argInts[i][2], argInts[i][3]);
//				System.out.println("\ttwoTailedP = " + twoTailedP);
//			}
//		}

	}


	static public class SmartFileReader extends FileReader {

		public SmartFileReader(File file) throws FileNotFoundException {
			super(file);
		}

		public SmartFileReader(FileDescriptor fd) {
			super(fd);
		}

		public SmartFileReader(String fileName) throws FileNotFoundException {
			super(fileName);
		}

		public boolean ready() throws IOException {
			super.ready();
			// Implemented because sun.nio.cs.StreamDecoder doesn't implement ready() properly.
			return true;
		}

	}

	final static double MIN_FREQ_FOR_HOM = 0.70;

	/**
	 * Runs the main execution logic
	 * @param args		Command-line arguments
	 */
	public static void main(String[] args) {

		String usage = "VarScan v2.4.4\n\n***NON-COMMERCIAL VERSION***\n\nUSAGE: java -jar VarScan.jar [COMMAND] [OPTIONS] \n\n";
		usage = usage + "COMMANDS:\n" +
				"\tcopynumber\t\t\tDetermine relative tumor copy number from tumor-normal pileups\n" +
				"\n";

		if(args.length > 0)
		{
			HashMap<String, String> params = getParams(args);

			if(args[0].equals("copynumber"))
			{
				copynumber(args, params);
			}

			else
			{
				System.err.println("Command not recognized\n" + usage);
			}
		}
		else
		{
			System.err.println(usage);
		}

	}


	/**
	 * Calls SNPs from a pileup file
	 *
	 * @param	args			Command-line arguments and parameters
	 * @param	callType		"SNP", "INDEL", or "CNS"
	 */


	/**
	 * Calls SNPs from an mpileup file
	 *
	 * @param	args			Command-line arguments and parameters
	 * @param	callType		"SNP", "INDEL", or "CNS"
	 */


	/**
	 * Obtains read counts for a list of variants
	 *
	 * @param	args	Command-line arguments
	 */



	/**
	 * Calls somatic/germline/LOH variants from normal and tumor pileup files
	 *
	 * @param	args	Command-line arguments
	 */


	/**
	 * Calls SNPs in a father-mother-child trio from an mpileup file
	 *
	 * @param	args			Command-line arguments and parameters
	 * @param	callType		"SNP", "INDEL", or "CNS"
	 */



	/**
	 * Determines tumor copy number from normal and tumor pileup files
	 *
	 * @param	args	Command-line arguments
	 */
	public static void copynumber(String[] args, HashMap<String, String> params)
	{
		if(params.containsKey("mpileup"))
		{
			 Copynumber myCopynumber = new Copynumber(args, true);
		}
		else
		{
			 Copynumber myCopynumber = new Copynumber(args);
		}

	}


	/**
	 * Filters variants by coverage, significance, frequency, etc.
	 *
	 * @param	args	Command-line arguments
	 */



	/**
	 * Applies false positive filter using bam-readcount information
	 *
	 * @param	args	Command-line arguments
	 */


	/**
	 * Filters variants by coverage, significance, frequency, etc.
	 *
	 * @param	args	Command-line arguments
	 */


	/**
	 * Splits VarScan output according to somatic status and confidence
	 *
	 * @param	args	Command-line arguments
	 */


	/**
	 * Calls somatic copy number events from copynumber output
	 *
	 * @param	args	Command-line arguments
	 */


	/**
	 * Compares two lists of positions/variants
	 *
	 * @param	args	Command-line arguments
	 */



	/**
	 * Limits pileup or variant files to a list of positions or regions
	 *
	 * @param	args	Command-line arguments
	 */

	/**
	 * Reports region coverage from a BAM file
	 *
	 * @param	args	Command-line arguments
	 */




	/**
	 * Parses and verifies any command-line parameters
	 *
	 * @param	args	Command-line arguments
	 * @return			HashMap of parameter names and their values
	 */
	static HashMap getParams(String[] args)
	{
		HashMap<String, String> params = new HashMap<String, String>();

		// Parse out command line arguments //

		String arg = "";
		String value = "";
		int i = 0, j = 0;

		// Go through each argument in the command line //

		while (i < args.length)
		{
			j = i + 1;
			arg = args[i];

			// If the argument starts with a hyphen, make use of it //

			if (arg.startsWith("-"))
			{
				// Remove leading hyphens //
				while(arg.startsWith("-"))
				{
					arg = arg.replaceFirst("-", "");
				}

				// Parse out parameters followed by values //

				if (i < args.length && j < args.length && !args[j].startsWith("-"))
				{
					value = args[j];
					params.put(arg, value);
				}

				// Set other parameters to true //

				else
				{
					params.put(arg, "true");
				}
			}

			i++;
		}

		return(params);
	}


	/**
	 * Gets the infile from command line or input buffer
	 *
	 * @param	args	Command-line arguments
	 * @return			HashMap of parameter names and their values
	 */
	static BufferedReader getInfile(String[] args)
	{
		BufferedReader in = null;

	    try
	    {
	    	// Declare file-parsing variables //

	    	String line;

	    	// Check for file on command line //

	    	if(args.length > 1 && !args[1].startsWith("-"))
	    	{
	    		File infile = new File(args[1]);
	    		if(infile.exists())
	    		{
	    			// Parse the infile //
	    			System.err.println("Reading input from " + args[1]);
	    			in = new BufferedReader(new SmartFileReader(args[1]));
	    		}
	    		else
	    		{
//    				System.err.println("File not found: " + args[1] + "\n");
//    				System.exit(10);
	    		}
	    	}

	    	// If no file from command line was parsed, try for piped input //

	    	if(in == null)
	    	{
		    	// Check the input stream //
		    	InputStreamReader instream = new InputStreamReader(System.in);
		    	Thread.sleep(1000);

		    	int num_naps = 0;

	    		while(!instream.ready())
	    		{
	    			System.err.println("Input stream not ready, waiting for 5 seconds...");
	    			Thread.sleep(5000);
	    			num_naps++;

	    			if(num_naps >= 100)
	    			{
	    				System.err.println("ERROR: Gave up waiting after 500 seconds...\n");
	    				System.exit(10);
	    			}
	    		}

		    	// If we have piped input, proceed with it //

		    	if(instream.ready())
		    	{
		    		System.err.println("Reading input from STDIN");
			    	in = new BufferedReader(instream);
		    	}
	    	}
	    }
	    catch(Exception e)
	    {
	    	System.err.println("ERROR: Unable to open input stream\n");
	    	System.exit(10);
	    }

		return(in);
	}


	/**
	 * Counts the depth of read bases meeting a minimum quality
	 *
//	 * @param	refBase		Reference base at this position
//	 * @param	readBases	String of read bases from pileup
	 * @param	readQuals	String of read base qualities from pileup
	 * @param	minAvgQual	Integer of minimum required base quality to count a base.
	 * @return	results		HashMap<String, String> of results for each allele
	 */
	static int qualityDepth(String readQuals, int minAvgQual)
	{
		int baseQuality = 0;
		int qualityDepth = 0;

		char[] arrQualities = readQuals.toCharArray();

		// Set quality position offset //
		int j = 0;

		// Go through each base //

		for(j = 0; j < arrQualities.length; j++)
		{
				baseQuality = arrQualities[j] - 33;
				if(baseQuality >= minAvgQual)
				{
					qualityDepth++;
				}
		}

		return(qualityDepth);
	}


	/**
	 * Calculates significance of read counts between two samples
	 *
	 * @param	expReads1	Reads supporting allele 1 (expected)
	 * @param	expReads2	Reads supporting allele 2 (expected)
	 * @param	obsReads1	Reads supporting allele 1 (observed)
	 * @param	obsReads2	Reads supporting allele 2 (observed)
	 * @return	p-value 	P-value from Fisher's Exact Test
	 */
	public static double getSignificance(int expReads1, int expReads2, int obsReads1, int obsReads2)
	{
		double pValue = 1;

		if(expReads1 < 0)
			expReads1 = 0;

		if(expReads2 < 0)
			expReads2 = 0;

		if(obsReads1 < 0)
			obsReads1 = 0;

		if(obsReads2 < 0)
			obsReads2 = 0;

		// Set up fisher's exact test //

		FishersExact fisher = new FishersExact(expReads1 + expReads2 + obsReads1 + obsReads2 + 100);

		// Calculate a p-value //

		pValue = fisher.getRightTailedP(expReads1, expReads2, obsReads1, obsReads2);
		int fisher_max = 1000;
		int num_tries = 0;

		while(Double.isNaN(pValue) && num_tries < 10)
		{
			fisher = new FishersExact(expReads1 + expReads2 + obsReads1 + obsReads2 + fisher_max);
			//pValue = fisher.getTwoTailedP(expReads1, expReads2, obsReads1, obsReads2);
			pValue = fisher.getRightTailedP(expReads1, expReads2, obsReads1, obsReads2);
			fisher_max = fisher_max + 1000;
			num_tries++;
		}

		if(num_tries >= 10)
			System.err.println("Warning: unable to calculate p-value failure: " + expReads1 + "," + expReads2 + "," + obsReads1 + "," + obsReads2);

		// If p-value is 1, do left-sided test //

		if(pValue >= 0.999)
		{
			pValue = fisher.getLeftTailedP(expReads1, expReads2, obsReads1, obsReads2);

			while(Double.isNaN(pValue))
			{
				fisher = new FishersExact(expReads1 + expReads2 + obsReads1 + obsReads2 + fisher_max);
				//pValue = fisher.getTwoTailedP(expReads1, expReads2, obsReads1, obsReads2);
				pValue = fisher.getLeftTailedP(expReads1, expReads2, obsReads1, obsReads2);
				fisher_max = fisher_max + 1000;
			}
		}

		return(pValue);
	}


}

package crispr;

import java.io.File;
import java.io.FileReader;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Locale;
import java.util.StringTokenizer;
import org.apache.commons.cli.*;


public class CrisprDist {
	static CommandLineParser parser = new GnuParser();
	static CommandLine cmd;
	static Options options = new Options();
	static HelpFormatter formatter = new HelpFormatter();
	

	static Option outputD   = OptionBuilder.withArgName( "directory" )
	.hasArg()
	.withDescription(  "output directory" )
	.create( "outputdir" );
	
	static Option cassettesD   = OptionBuilder.withArgName( "directory" )
	.hasArg()
	.withDescription(  "directory with cassettes" )
	.create( "cassettes" );
	static Option names_table   = OptionBuilder.withArgName( "file" )
		.hasArg()
		.withDescription(  "file with table of alternative names(see format below). if 'numbers' exist then creates file with such name" )
		.create( "names_table" );
		
	static Option gbopt   = OptionBuilder.withArgName( "float" )
	.hasArg()
	.withType("float")
	.withDescription(  "begining gap score (default -3)" )
	.create( "gb" );
	
	static Option goopt   = OptionBuilder.withArgName( "float" )
	.hasArg()
	.withType("float")
	.withDescription(  "opening gap score (default -5)" )
	.create( "go" );
	
	static Option geopt   = OptionBuilder.withArgName( "float" )
	.hasArg()
	.withType("float")
	.withDescription(  "elongation gap score (default -1)" )
	.create( "ge" );
	
	static Option matchopt   = OptionBuilder.withArgName( "float" )
	.hasArg()
	.withType("float")
	.withDescription(  "match score (default 5)" )
	.create( "match" );
	
	
	
	//Weights:
	public static float gb = -3;
	public static float go = -5;
	public static float ge = -1;
	public static float match = 5;
	public static float mismatch = Float.NEGATIVE_INFINITY;
	
	//Main Objects to work with
	private static HashMap<String, ArrayList<Integer>> cassettes = new HashMap<String, ArrayList<Integer>>();
	private static HashMap<String, String> altNames = new HashMap<String, String>();
	private static ArrayList<String> spacersSet = new ArrayList<String>();
	private static ArrayList<String> names = new ArrayList<String>();
	
	//Order of matrices: 0 F, 1 H, 2 V
	//Matrices:
	private static float[][] F;
	private static float[][] H;
	private static float[][] V;
	private static int[][] piF;
	private static int[][] piH;
	private static int[][] piV;
	
	//Files and Directory
	private static File cassettesDir;
	private static File outputDir;
	private static PrintWriter namesTableFile;
	private static PrintWriter distancesFile;
	private static String namesTable= "NamesTable.txt";
	
	public static boolean inverse = false;
	
	//Options
	public static boolean outputNames = false;
	
	public static void main(String[] arg) throws Exception
	{
		
		
		options.addOption("help", false, "print this help");
		options.addOption(outputD);
		options.addOption(cassettesD);
		options.addOption(names_table);
		options.addOption(gbopt);
		options.addOption(goopt);
		options.addOption(geopt);
		options.addOption(matchopt);
		options.addOption("inverse", false, "inverts order of spacers(read below)");
		options.addOption("numbers", false, "make a table with alternative names as numbers");
		
		
	    try {
	        // parse the command line arguments
	        cmd = parser.parse( options, arg );
	    }
	    catch( ParseException exp ) {
	        // oops, something went wrong
	        System.err.println( "Parsing failed.  Reason: " + exp.getMessage() );
	    }
	    
	    	    
	    if(cmd.hasOption("help") || cmd.getOptions().length == 0){
	    	String helpheader = "This program estimates matrix of distances between CRISPR-cassettes\n";
	    	String helpfooter = "\n\tOrder:\n" +
	    			"Straight order means that the first spacer in cassette is the last added spacer" +
	    			"\n\n\t Alternative names file format:\n" +
	    			"File contains two colomns separated by tabulator. The first colomn corresponds " +
	    			"to names of files with cassettes. And the second one corresponds to altenative " +
	    			"names in output file.";
	    	formatter.printHelp( "java -jar crisprdist.jar", helpheader,options,helpfooter, true ); 
	    	return;
	    }
		
		//BufferedReader in = new BufferedReader(new InputStreamReader(System.in));
		try{
			cassettesDir = new File(cmd.getOptionValue("cassettes"));
			if(!cassettesDir.isDirectory()) throw new IllegalArgumentException("Destination directory ["
					+ cassettesDir.getPath() + "] is not valid.");
			
			/*System.out.print("Input weights: Begining Gap, Opening Gap, Elongation Gap, Match (default: -3 -5 -1 5): ");
			StringTokenizer st;
			while(true){
				st = new StringTokenizer(in.readLine(), " \t,");
				if(st.countTokens() == 0) {
					break; 
				}
				if(st.countTokens() == 4) {
					gb = Float.parseFloat(st.nextToken());
					go = Float.parseFloat(st.nextToken());
					ge = Float.parseFloat(st.nextToken());
					match = Float.parseFloat(st.nextToken());
					break; 
				}
				System.out.println("4 numbers are necessary, try again...");
			}*/
			if(cmd.hasOption("gb")) gb = Float.parseFloat(cmd.getOptionValue("gb"));
			if(cmd.hasOption("go"))go = Float.parseFloat(cmd.getOptionValue("go"));
			if(cmd.hasOption("ge"))ge = Float.parseFloat(cmd.getOptionValue("ge"));
			if(cmd.hasOption("match"))match = Float.parseFloat(cmd.getOptionValue("match"));
			
			//System.out.print("Direction: is order from Leader-sequence? (empty for 'yes') ");
			//if(!in.readLine().isEmpty()) inverse = true;
			if(cmd.hasOption("inverse")) inverse = true;
			for(File f: cassettesDir.listFiles()){
				BufferedReader br = new BufferedReader(new FileReader(f));
				cassettes.put(f.getName(), new ArrayList<Integer>());
				names.add(f.getName());
				String line = br.readLine();
				String spacer = "";
				int index;
				while(line != null){
					if(line.charAt(0) == '>'){
						spacer = "";
						while((line = br.readLine()) != null && line.charAt(0) != '>'){
							spacer += line.replaceAll("\\s", "");
						}
						index = spacersSet.indexOf(spacer);
						if(index == -1){
							cassettes.get(f.getName()).add(spacersSet.size());
							spacersSet.add(spacer);
						}else{
							cassettes.get(f.getName()).add(index);
						}
					}
				}
				br.close();	
				//System.out.println(cassettes.get(f.getName()).get(0)+"\t"+cassettes.get(f.getName()).get(cassettes.get(f.getName()).size()-1));
				
				
				if(inverse){
					int temp;
					for(int smaller = 0, bigger=cassettes.get(f.getName()).size() - 1; smaller < bigger; smaller++, bigger--){
						temp = cassettes.get(f.getName()).get(smaller);
						cassettes.get(f.getName()).set(smaller, cassettes.get(f.getName()).get(bigger));
						cassettes.get(f.getName()).set(bigger, temp);
					}
				}
				//System.out.println(cassettes.get(f.getName()).get(0)+"\t"+cassettes.get(f.getName()).get(cassettes.get(f.getName()).size()-1));
			}
		}catch(IllegalArgumentException e){
			System.out.println(e.getLocalizedMessage());
			return;
		}catch(Exception e){
			throw e;
		}	
		
		try {
			outputDir = cmd.hasOption("outputdir") ? new File(cmd.getOptionValue("outputdir")): new File(".");
			outputDir.mkdir();
			File distFile = new File(outputDir + "/" + cassettesDir.getName() + ".dist");
			if(distFile.exists()) distFile.delete();
			distancesFile = new PrintWriter(outputDir + "/" + cassettesDir.getName() + ".dist");
			distancesFile.print("\t" + names.size() + "\n");
			
			//Some options
			//System.out.print("Using alternative short names in format [Name\tAltName] from file (empty if no such file): ");
			//String userTable = in.readLine();
			//String numbersAsNames = "";
			if (cmd.hasOption("names_table")) {
				namesTable = cmd.getOptionValue("names_table");
				if(cmd.hasOption("numbers")){
					//System.out.print("Create a file [" + cassettesDir.getName() + namesTable +"] with numbers for AltName (empty for 'yes')? ");
					//numbersAsNames = in.readLine();
					namesTableFile = new PrintWriter(outputDir + "/" + namesTable);
				}else{
					//File userTableFile = new File(cmd.getOptionValue("names_table"));
					BufferedReader userTableReader = new BufferedReader(new FileReader(cmd.getOptionValue("names_table")));
					String line;
					StringTokenizer altNamesST;
					while((line = userTableReader.readLine()) != null){
						altNamesST = new StringTokenizer(line, "\t");
						if(altNamesST.countTokens() != 2) throw new IllegalArgumentException("You must input an AltNames file in format: [Name\tAltName].");
						altNames.put(altNamesST.nextToken(), altNamesST.nextToken());						
					}
				}
			}
			//boolean printNamesToFile = userTable.isEmpty() && numbersAsNames.isEmpty();
						
			int nameIterator = 0;
			String name;
			
			//Calculating and making a file
			for(int row = 0; row < names.size(); row++){
				nameIterator++;
				
				if(cmd.hasOption("names_table")){
					if(cmd.hasOption("numbers")){
				//if(userTable.isEmpty()){
					//if(numbersAsNames.isEmpty()){
						name = String.valueOf(nameIterator);
					}else{
						System.out.println(altNames);
						name = altNames.get(names.get(row)).replaceAll("[\\s:;,\\(\\)\\[\\]]", "+");
					}
				}else{
					name = names.get(row).replaceAll("[\\s:;,\\(\\)\\[\\]]", "+");
				}
				//distancesFile.printf("%-10.10s", nameIterator + (outputNames ? "|" + names.get(row).replaceAll("[\\s:;,\\(\\)\\[\\]]", "_") : ""));
				distancesFile.printf("%-10s", name);    //prints the whole name of any length
				if(cmd.hasOption("names_table") && cmd.hasOption("numbers")) namesTableFile.println(names.get(row) + "\t" + nameIterator);
				//System.out.println(nameIterator + " : " + names.get(row));
				
				int newline = 0;
				for(int col = 0; col < row; col++){
					Integer[] A = cassettes.get(names.get(row)).toArray(new Integer[cassettes.get(names.get(row)).size()]);
					Integer[] B = cassettes.get(names.get(col)).toArray(new Integer[cassettes.get(names.get(col)).size()]);
					calculateMatrices(A, B);
					int[] gaps = getGaps();
					
					float minWeight = (A.length + B.length) * gb;
					float maxWeight = max(new float[] {A.length, B.length})[0] * match;
					float dist = 1 - ((gaps[0]*gb + gaps[1]*(go + ge) + gaps[2]*ge + gaps[3]*match) - minWeight) / (maxWeight - minWeight);
					distancesFile.printf(Locale.US, "%10.6f", dist);
					
					//Check
					/*System.out.printf(Locale.US, "%s - %s, gb:%d, go:%d, ge:%d, mat:%d\n", names.get(row), names.get(col), gaps[0], gaps[1], gaps[2], gaps[3]);
					for (int i = 0; i < F.length; i++){
						for (int j = 0; j < F[0].length; j++)
						{					
							System.out.printf(Locale.US, "%10.2f", F[i][j]);
						}
						System.out.println("");
					}
					System.out.println("");
					for (int i = 0; i < F.length; i++){
						for (int j = 0; j < F[0].length; j++)
						{					
							System.out.printf(Locale.US, "%10.2f", H[i][j]);
						}
						System.out.println("");
					}
					System.out.println("");
					for (int i = 0; i < F.length; i++){
						for (int j = 0; j < F[0].length; j++)
						{					
							System.out.printf(Locale.US, "%10.2f", V[i][j]);
						}
						System.out.println("");
					}*/
					newline++;
					if ((newline % 6) == 0) distancesFile.print("\n");
				}
				for(int col = row; col < names.size(); col++){
					distancesFile.printf(Locale.US, "%10.1f", 0.0);
					
					newline++;
					if ((newline % 6) == 0) distancesFile.print("\n");
				}
				
				distancesFile.print("\n");
			}
			distancesFile.close();
			if(cmd.hasOption("names_table") && cmd.hasOption("numbers")) namesTableFile.close();
		} catch (IOException e) {
			throw e;
		}
			
	}
	
	public static int[] getGaps()
	{
		int[] numOfGaps = {0, 0, 0, 0};
		
		int i = F.length - 1;
		int j = F[0].length - 1;
		int matrixNum = (int)max(new float[] {F[i][j], H[i][j], V[i][j]})[1];
		
		all:
		while(i >= 0 && j >= 0){
			switch(matrixNum){
			case 0:
				if(piF[i][j] == 3){
					numOfGaps[0] = i + j;
					break all;
				}
				numOfGaps[3]++;
				matrixNum = piF[i][j];
				i--;
				j--;
				break;
			case 1:
				if(piH[i][j] == 1){
					numOfGaps[2]++;
				}else{
					numOfGaps[1]++;
				}
				matrixNum = piH[i][j];
				i--;
				break;
			case 2:
				if(piV[i][j] == 2){
					numOfGaps[2]++;
				}else{
					numOfGaps[1]++;
				}
				matrixNum = piV[i][j];
				j--;
			}
		}
		return numOfGaps;
	}
	
	private static void calculateMatrices(Integer[] source, Integer[] dest)
	{
		F = new float[source.length+1][dest.length+1];
		V = new float[source.length+1][dest.length+1];
		H = new float[source.length+1][dest.length+1];
		piF = new int[source.length+1][dest.length+1];
		piV = new int[source.length+1][dest.length+1];
		piH = new int[source.length+1][dest.length+1];
		
				
		F[0][0] = 0;
		piF[0][0] = 3;
		H[0][0] = V[0][0] = Float.NEGATIVE_INFINITY;
		for (int i = 1; i < source.length + 1; i++){
			H[i][0] = V[i][0] = Float.NEGATIVE_INFINITY;
			F[i][0] = i*gb;
			piF[i][0] = 3;
		}
		for (int j = 1; j < dest.length + 1; j++){
			H[0][j] = V[0][j] = Float.NEGATIVE_INFINITY;
			F[0][j] = j*gb;
			piF[0][j] = 3;
		}
		
		float ff, fh, fv, init; 		//Values for F
		float hf, hv, hh; 				//Values for H
		float vf, vh, vv;				//Values for V
		
		float[] max;
		
		for (int i = 1; i < source.length + 1; i++)
			for (int j = 1; j < dest.length + 1; j++)
				{					
					ff = F[i-1][j-1] + simular(source[i-1] , dest[j-1]);
					fh = H[i-1][j-1] + simular(source[i-1] , dest[j-1]);
					fv = V[i-1][j-1] + simular(source[i-1] , dest[j-1]);
					init = gb * (i + j);
					max = max(new float[] {ff, fh, fv, init});
					F[i][j] = max[0];
					piF[i][j] = (int)max[1];
					
					hf = piF[i-1][j] != 3 ?  F[i-1][j] + go + ge : Float.NEGATIVE_INFINITY;
					hv = V[i-1][j] + go + ge;
					hh = H[i-1][j] + ge;
					max = max(new float[] {hf, hh, hv});
					H[i][j] = max[0];
					piH[i][j] = (int)max[1];
					
					vf = piF[i][j-1] != 3 ? F[i][j-1] + go + ge : Float.NEGATIVE_INFINITY;
					vh = H[i][j-1] + go + ge;
					vv = V[i][j-1] + ge;
					max = max(new float[] {vf, vh, vv});
					V[i][j] = max[0];
					piV[i][j] = (int)max[1];
				}
	}

	private static float simular(int first, int second)
	{
		return (first == second ? match : mismatch);
	}
	
	private static float[] max(float[] args){
		float[] result = {args[0], 0};
		for(int i = 1; i < args.length; i++){
			if(args[i] > result[0]){
				result[0] = args[i];
				result[1] = i;
			}
		}
		return result;
	}
			
}

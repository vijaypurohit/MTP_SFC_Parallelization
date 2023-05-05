/*!
 * SND Library Network Reader.
 * Original Author: Prateekshya Priyadarshini (Thanks to her)
 * Modified: Vijay Purohit (for requirements according to my code)
*/
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.InputStreamReader;
import java.io.IOException;

import java.io.FileReader;
import java.io.FileWriter;
import java.io.File;

import java.util.List;
import java.util.ArrayList;
import java.util.Map;
import java.util.LinkedHashMap;
import java.util.TreeMap;
import java.util.Set;
import java.util.LinkedHashSet;
import java.util.Arrays;

class ImportNetwork
{
    public final static double AVERAGE_RADIUS_OF_EARTH_KM = 6371;
    public static double calculateDistanceInKilometer(double lat1, double lng1, double lat2, double lng2) {

        double latDistance = Math.toRadians(lat1 - lat2);
        double lngDistance = Math.toRadians(lng1 - lng2);

        double a = Math.sin(latDistance / 2) * Math.sin(latDistance / 2)
          + Math.cos(Math.toRadians(lat1)) * Math.cos(Math.toRadians(lat2))
          * Math.sin(lngDistance / 2) * Math.sin(lngDistance / 2);

        double c = 2 * Math.atan2(Math.sqrt(a), Math.sqrt(1 - a));

        return (AVERAGE_RADIUS_OF_EARTH_KM * c);
    }//calculateDistanceInKilometer

	public static void main(String[] args) throws IOException
	{
        BufferedReader reader = new BufferedReader( new InputStreamReader(System.in));
        String filepath;
        System.out.println("Enter SNDLib filepath /files_input/ ");
        filepath = reader.readLine(); // taking string input
        String[] fileTokens = filepath.trim().split("/+");

		boolean nodeReading = false, linkReading = false;
		Map<String,List<Double>> NODES = new LinkedHashMap<>();
		Set<String> LINKS = new LinkedHashSet<>();
		BufferedReader fileReader = new BufferedReader(new FileReader(new File(filepath))); //files_input/NewYork/
		String line = "";
		while ((line = fileReader.readLine()) != null)
		{
			if (line.equals("NODES (")) {
				nodeReading = true;
				continue;
			}
			else if (line.equals("LINKS (")) {
				linkReading = true;
				continue;
			}
			else if (line.equals(")")){
				nodeReading = linkReading = false;
			}
			else if (line.equals("# DEMAND SECTION")){
				break;
			}

			if (nodeReading) {
				String[] tokens = line.trim().split(" +");
				NODES.put(tokens[0],Arrays.asList((double)(NODES.size()+1),Double.parseDouble(tokens[2]),Double.parseDouble(tokens[3])));
			}

			if (linkReading) {
				String[] tokens = line.trim().split(" +");
				LINKS.add(tokens[2]+" "+tokens[3]);
			}
		} 
		fileReader.close();

		int numV = NODES.size(); ///< number of vertices in the network
		int numE = LINKS.size(); ///< number of edges in the network.
		
		int[] degrees = new int[numV+1];
		double[][] network = new double[numV+1][numV+1]; //< adj matrix
		StringBuilder l = new StringBuilder();
		int lsize=1;
		for (String link : LINKS)
		{
			String[] tokens = link.split(" ");
			List<Double> node1 = NODES.get(tokens[0]); // src node
			List<Double> node2 = NODES.get(tokens[1]); // dst node
			int n1 = node1.get(0).intValue(); // src index value
			int n2 = node2.get(0).intValue(); // dst index value
			double x1 = node1.get(1); // lat of src
			double x2 = node2.get(1); // lat of dst
			double y1 = node1.get(2); // long of src
			double y2 = node2.get(2); // long of dst
			network[n1][n2] = network[n2][n1] = calculateDistanceInKilometer(x1,y1,x2,y2)*1000 ; //Math.sqrt(Math.pow((x2-x1),2)+Math.pow((y2-y1),2));
			l.append(n1 + " "+n2 +" "+ network[n1][n2] );
			if(lsize != LINKS.size()) l.append("\n");
			lsize++;
			++degrees[n2];
			++degrees[n1];
		}

		BufferedWriter fileWriter = new BufferedWriter(new FileWriter(new File(fileTokens[0]+"/network_"+fileTokens[0]+".txt"))); //files_input/NewYork/
		// fileWriter.write("Network Data");fileWriter.newLine();
		
		fileWriter.write(numV + " L\n");
		fileWriter.write(numV + " \n"); 
		for(Map.Entry node: NODES.entrySet()){

				int nodeIndex = ((Double)((List)node.getValue()).get(0)).intValue();
				int numCores = (degrees[nodeIndex]);//(degrees[nodeIndex]>=8 ? 4: 2);
 
				fileWriter.write(nodeIndex + " " + numCores);fileWriter.newLine();
		}		
		

		fileWriter.write(l.toString());

		fileWriter.close();
        System.out.println(filepath+" converted.");
	}//main
}//ImportNetwork Class
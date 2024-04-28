package progistar.scan.data;

import java.util.ArrayList;
import java.util.Hashtable;

public class LocTable {
	public Hashtable<String, Hashtable<String, LocationInformation>> table = new Hashtable<String, Hashtable<String, LocationInformation>>(); 

	public boolean putLocation (LocationInformation lInfo) {
		Hashtable<String, LocationInformation> gTable = table.get(lInfo.inputSequence);
		if(gTable == null) {
			gTable = new Hashtable<String, LocationInformation>();
			table.put(lInfo.inputSequence, gTable);
		}
		
		if(gTable.get(lInfo.getKey()) == null) {
			gTable.put(lInfo.getKey(), lInfo);
			return true;
		}
		
		return false;
	}
	
	public ArrayList<String> getLocations (String inputSequence) {
		ArrayList<String> locations = new ArrayList<String>();
		String na = "Not found";
		
		Hashtable<String, LocationInformation> gTable = table.get(inputSequence);
		if(gTable != null) {
			gTable.forEach((key, info)->{
				locations.add(info.getRes());
			});
		}
		
		if(locations.size() == 0) {
			locations.add(na);
		}
		
		return locations;
	}
}

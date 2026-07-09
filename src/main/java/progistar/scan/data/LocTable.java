package progistar.scan.data;

import java.util.ArrayList;
import java.util.HashMap;

public class LocTable {
	public HashMap<String, HashMap<String, LocationInformation>> table = new HashMap<String, HashMap<String, LocationInformation>>(); 

	public boolean putLocation (LocationInformation lInfo) {
		HashMap<String, LocationInformation> gTable = table.get(lInfo.inputSequence);
		if(gTable == null) {
			gTable = new HashMap<String, LocationInformation>();
			table.put(lInfo.inputSequence, gTable);
		}
		
		LocationInformation slInfo = gTable.get(lInfo.getKey());
		if(slInfo == null) {
			gTable.put(lInfo.getKey(), lInfo);
			return true;
		} else {
			lInfo.readCounts.forEach((barcodeId, value)->{
				Long val = slInfo.readCounts.get(barcodeId);
				if(val == null) {
					val = 0L;
				}
				slInfo.readCounts.put(barcodeId, (value + val));
			});
		}
		
		return false;
	}
	
	public ArrayList<LocationInformation> getLocations (String inputSequence) {
		ArrayList<LocationInformation> locations = new ArrayList<LocationInformation>();
		HashMap<String, LocationInformation> gTable = table.get(inputSequence);
		if(gTable != null) {
			gTable.forEach((key, info)->{
				locations.add(info);
			});
		}
		
		// if there is no matched reads
		if(locations.size() == 0) {
			LocationInformation nullLocation = new LocationInformation();
			nullLocation.readCounts.put(Constants.DEFAULT_BARCODE_ID, 0L);
			nullLocation.inputSequence = inputSequence;
			nullLocation.obsPeptide = inputSequence;
			locations.add(nullLocation);
		}
		
		return locations;
	}
}

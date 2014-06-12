package assignment4.evolution;

import java.util.Random;

/**
 * {@code MastDistribution}s are distributions of {@code TransmitterMast}s
 * which are stored as one map of masts and one of resulting signal strengths.
 * This class offers functions for placing and removing masts and tells the
 * total cost and the number of sufficiently covered areas.</br>
 * 
 * Multiple {@code MastDistribution}s form a {@code Generation}.</br></br>
 * 
 * The goal is to find a distribution that provides at least 90% of all areas
 * with a signal strength of 3 or above and is as inexpensive as possible.
 * 
 * @author Andreas Burmeister
 * 
 * @version 0.3.2 05/27/14
 * 
 * @see TransmitterMast
 * @see Generation
 */

public class MastDistribution {
	
	private long overallCost;
	//The total price for the configuration.
	
	public long getOverallCost() {
		return this.overallCost;
	}
	
	private int numberOfMasts;
	//The number of masts currently used.
	
	public int getNumberOfMasts() {
		return this.numberOfMasts;
	}
	
	private TransmitterMast[][] mastMap =
		new TransmitterMast[CountryMap.mapHeight][CountryMap.mapWidth];
	//A two-dimensional array to store where which masts are placed.
	
	public TransmitterMast[][] getMastMap() {
		return this.mastMap;
	}
	
	private int[][] coverageMap =
		new int[CountryMap.mapHeight][CountryMap.mapWidth];
	//A two-dimensional array to store how well which areas are covered.
	
	public int[][] getCoverageMap() {
		return this.coverageMap;
	}
	
	public int coverage() {
	//Tells how many areas are covered with sufficient signal strength (>=3).
		int coverage = 0;
		int[][] coveragemap = this.coverageMap;
		for(int row = 0; row < CountryMap.mapHeight; row++) {
			for(int col = 0; col < CountryMap.mapWidth; col++) {
				if (CountryMap.countryMap[row][col] == 1)
					coverage += (coveragemap[row][col] >= 3 ? 1 : 0);
			}
		}
		return coverage;
	}
	
	public void addMast(TransmitterMast mast, int posX, int posY)
	throws Exception {
	//Places a mast at the given coordinates.
		if (!(posY<0||posY>=CountryMap.mapHeight||
			posX<0||posX>=CountryMap.mapWidth)) {
				if (!(CountryMap.countryMap[posY][posX] == 0)) {
					if (this.mastMap[posY][posX] == null) {
						this.mastMap[posY][posX] = mast;
						int[][] coverage = mast.getAreaOfCoverage();
						for (int row = 0; row < 9; row++) {
							for (int col = 0; col < 9; col++) {
								if (!(posY+row-4<0||
									posY+row-4>=CountryMap.mapHeight||
									posX-col+4<0||
									posX-col+4>=CountryMap.mapWidth))
										this.coverageMap
											[posY+row-4][posX-col+4]
											+= coverage[row][col];
							}
						}
						this.overallCost+=mast.getCost();
						this.numberOfMasts++;
					} else {
						throw new Exception(
							"Cannot add another mast in this place."+
							"(X="+posX+", Y="+posY+")");
					}
				} else {
					throw new Exception(
						"Cannot add a mast outside the country."+
						"(X="+posX+", Y="+posY+")");
				}
		} else {
			throw new Exception(
				"Cannot add a mast outside the map."+
				"(X="+posX+", Y="+posY+")");
		}
	}
	
	public void removeMast(int posX, int posY)
	throws Exception {
	//Removes a mast from the given coordinates.
		if (!(posY<0||posY>=CountryMap.mapHeight||
			posX<0||posX>=CountryMap.mapWidth)) {
				if (!(this.mastMap[posY][posX] == null)) {
					this.overallCost-=mastMap[posY][posX].getCost();
					int[][] coverage = mastMap[posY][posX].getAreaOfCoverage();
					for (int row = 0; row < 9; row++) {
						for (int col = 0; col < 9; col++) {
							if (!(posY+row-4<0||
								posY+row-4>=CountryMap.mapHeight||
								posX-col+4<0||
								posX-col+4>=CountryMap.mapWidth))
									this.coverageMap[posY+row-4][posX-col+4] -=
										coverage[row][col];
						}
					}
					this.mastMap[posY][posX] = null;
					this.numberOfMasts--;
				} else {
					throw new Exception(
						"Cannot remove a non-existent mast."+
						"(X="+posX+", Y="+posY+")");
				}
		} else {
			throw new Exception(
				"Cannot remove a mast outside the map."+
				"(X="+posX+", Y="+posY+")");
		}
	}
	
	public TransmitterMast.MastType checkMast(int posX, int posY)
	throws Exception {
	//Tells what type of mast is placed at the given coordinates.
		if (!(posY<0||posY>=CountryMap.mapHeight||
			posX<0||posX>=CountryMap.mapWidth)) {
				if (!(this.mastMap[posY][posX] == null)) {
					return this.mastMap[posY][posX].getMastType();
				} else {
					throw new Exception(
						"Cannot check for a non-existent mast."+
						"(X="+posX+", Y="+posY+")");
				}
		} else {
			throw new Exception(
				"Cannot check for a mast outside the map."+
				"(X="+posX+", Y="+posY+")");
		}
	}
	
	public void clear() {
	//Deletes all entries of masts for this MastDistribution.
		for (int row = 0; row < CountryMap.mapHeight; row++) {
			for (int col = 0; col < CountryMap.mapWidth; col++) {
				this.mastMap[row][col] = null;
				this.coverageMap[row][col] = 0;
				this.overallCost = 0;
				this.numberOfMasts = 0;
			}
		}
	}
	
	public void randomize(int n) {
	//Fills this MastDistribution randomly with up to TransmitterMasts.
		TransmitterMast tma = new TransmitterMast(TransmitterMast.MastType.A);
	    TransmitterMast tmb = new TransmitterMast(TransmitterMast.MastType.B);
	    TransmitterMast tmc = new TransmitterMast(TransmitterMast.MastType.C);
	    this.clear();

	    Random rnd = new Random();
		int masts = rnd.nextInt(100)+1;
		for (int i = 0; i < masts; i++) {
			TransmitterMast mast = tma;
			switch (rnd.nextInt(3)) {
				case 0:
					mast = tma;
					break;
				case 1:
					mast = tmb;
					break;
				default:
					mast = tmc;
					break;
			}
			try {
				this.addMast(mast,
					(int)(Math.random()*CountryMap.mapWidth-1),
					(int)(Math.random()*CountryMap.mapHeight-1));
			} catch (Exception ex) {
			}
		}
	}
}
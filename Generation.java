package assignment4.evolution;

/**
 * A {@code Generation} consists of a number of individuals, in this case
 * {@code MastDistribution}s. This class offers some basic functionality such
 * as randomization to work with a Generation (as in {@code Evolution}).
 * 
 * @author Andreas Burmeister
 * 
 * @version 0.2.2 05/27/14
 * 
 * @see MastDistribution
 * @see Evolution
 */

public class Generation {
	
	private MastDistribution[] individuals;
	//An array of individual MastDistributions forms the generation.
	
	public MastDistribution[] getIndividuals() {
		return this.individuals;
	}
	
	public void setIndividual
		(MastDistribution distribution, int individualNr) {
	//Replaces an individual of the generation with another one.
		this.individuals[individualNr] = distribution;
	}
	
	public void randomize(int individualNr, int n) {
	//Randomizes one individual of this generation.
		this.individuals[individualNr] = new MastDistribution();
		this.individuals[individualNr].randomize(n);
	}
	
	public void randomizeAll(int n) {
	//Randomizes all individuals of this generation.
		for (int i = 0; i < this.individuals.length; i++) {
			this.individuals[i] = new MastDistribution();
			this.individuals[i].randomize(n);
		}
	}
	
	public Generation(int size) {
		this.individuals = new MastDistribution[size];
		for (int i = 0; i < size; i++) {
			this.individuals[i] = new MastDistribution();
		}
	}

	public int size() {
		return this.individuals.length;
	}
}
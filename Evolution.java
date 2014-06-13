package assignment4.evolution;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Random;
import java.util.Set;



/**
 * This is the class with the {@code main} method - it initializes the first
 * generation and then evolves it iteratively.</br>
 * To run the program, follow these steps:
 * 	Run -> Run Configurations... -> Java Application -> Right Click
 * 	-> New -> Main Class: assignment4.evolution.Evolution -> Run</br></br>
 * 
 * Your task is to implement the method {@code evolve()}.</br>
 * It has to manage the transition from one generation to the next.</br>
 * Therefore you will first need to rate the individuals (MastDistibutions)
 * based on their fitness (how well they fulfill the requirements).</br>
 * Second, you have to create a new generation from the old one by using:
 * <ul>
 * <li>selection: fit individuals survive and make it to the next generation
 * <li>crossover: features of fit individuals are recombined
 * <li>mutation: random changes occur, thus varying the genetic material
 * </ul>
 * 
 * This literature might help you get an overview of genetic algorithms or
 * evolutionary programming in general:</br>
 * 
 * <ul>
 * <li>Teske 2013: Computerprogramme züchten. In HPImgzn issue 14, pp. 50-53
 * <li>Saake, Sattler: Algorithmen und Datenstrukturen, 5th edition, pp. 86-88
 * <li>{@link http://en.wikipedia.org/wiki/Genetic_algorithm}
 * </ul>
 *
 * @author Andreas Burmeister
 * 
 * @version 0.3 05/27/14
 * 
 * @see Generation
 * @see MastDistribution
 * @see TransmitterMast
 */

public class Evolution {
	private final long budget;
	private final float minCoverage;
	private MastDistribution solution;
	private Generation currentGeneration;
	private int generationNumber;

	public Evolution(int generationSize, long budget, float minCoverage) {
		currentGeneration = new Generation(generationSize);
	    currentGeneration.randomizeAll(100); //default: 100
	    generationNumber = 1;
	    this.budget = budget;
	    this.minCoverage = minCoverage;
	    check();
	}
	
	public Generation getCurrentGeneration() {
		return currentGeneration;
	}

	public MastDistribution getSolution() {
		return solution;
	}

	public int getGenerationNumber() {
		return generationNumber;
	}
	
	private void check() {
		MastDistribution[] individuals = currentGeneration.getIndividuals();
		for (int i = 0; i < individuals.length; i++) {
			if (individuals[i].coverage() > minCoverage*CountryMap.countrySize &&
				individuals[i].getOverallCost() < budget) {
				solution = individuals[i];
				int a = solution.coverage();
				int b = CountryMap.countrySize;
				System.out.println(
					"HEUREKA!\n" +
					"Generation: " + generationNumber +
					", Individual: " + (i+1) + "\n" +
					"Coverage(3+): " + a + "/" + b + " (" +
					(int)(((float)a/b)*100) + " %)\n" +
					"Cost: " + solution.getOverallCost()
				);
				break;
			}
		}
	}
	
	public void evolve() {
		evolve(true);
	}
	
	public void evolve(boolean output) {
		///////////////////////////////////////////////////////////////////////
		//
		// TODO: implement this method (and maybe helper methods)
		//
		// How does the transition from one generation to the next happen?
		// Which individuals qualify to survive? How do you rate individuals?
		// Will mutations occur? Which kinds? Will there be crossovers? 
		//
		// e.g.:
		// crossover(currentGeneration.getIndividuals()[i], currentGeneration.getIndividuals()[j]);
		// mutate(currentGeneration.getIndividuals()[i]);
		// rateFitness(currentGeneration.getIndividuals()[i]);
		// 
		// Why do you choose to do what you do instead of something else?
		// Please comment everything you do. Help us understand your ideas.
		//
		///////////////////////////////////////////////////////////////////////
		
		// Get current and create next set of individuals
		ArrayList<MastDistribution> current = new ArrayList<MastDistribution>(Arrays.asList(currentGeneration.getIndividuals()));
		ArrayList<MastDistribution> next = new ArrayList<MastDistribution>();
		
		// Sort current generation by fitness
		Collections.sort(current, new Comparator<MastDistribution>() {
			public int compare(MastDistribution lhs, MastDistribution rhs) {
				return new Float(rateFitness(rhs)).compareTo(new Float(rateFitness(lhs)));
			}
		});
				
		// Keep some of the best
		for (int i = 0; i < current.size() / 6; ++i)
			next.add(current.get(i));
		
		// Mutate multiple times
		for (int i = 0; i < current.size(); ++i)
			for (int j = 0; j < current.size(); ++j)
				next.add(mutate(current.get(i), generationNumber));
		
		// Cross over with each other
		for (int i = 0; i < current.size(); ++i)
			for (int j = 0; j < current.size(); ++j)
				next.add(crossover(current.get(i), current.get(j)));
		
		// Cross over with the best
		for (int i = 0; i < current.size(); ++i)
			next.add(crossover(current.get(i), current.get(0)));
		
		// Cross over with random
		for (int i = 0; i < current.size(); ++i) {
			MastDistribution individual = new MastDistribution();
			individual.randomize(100);
			next.add(crossover(current.get(i), individual));
		}
		
		// Add some random
		for (int i = 0; i < current.size(); ++i) {
			MastDistribution individual = new MastDistribution();
			individual.randomize(100);
			next.add(individual);
		}
		
		// Sort candidates for next generation by fitness
		Collections.sort(next, new Comparator<MastDistribution>() {
			public int compare(MastDistribution lhs, MastDistribution rhs) {
				return new Float(rateFitness(rhs)).compareTo(new Float(rateFitness(lhs)));
			}
		});

		// Print stats of best individual
		if (output) {
			MastDistribution best = next.get(0);
			float coverage = (float) best.coverage() / CountryMap.countrySize;
			float cost = (float) best.getOverallCost() / budget;
			float rate = 0.31618f / (float) Math.pow(1 + generationNumber / 1000.0, 2.0);
			System.out.println("Generation: " + generationNumber + ", Coverage: " + coverage +  ", Cost: " + cost + ", Rate: " + rate);
		}
		
		// Create next generation form same amount of best individuals
		Generation nextGeneration = new Generation(current.size());
		for (int i = 0; i < current.size() && i < next.size(); ++i)
			nextGeneration.setIndividual(next.get(i), i);
		currentGeneration = nextGeneration;
		generationNumber++;
		if (output)
			check();
	}
	
	public float rateFitness(MastDistribution individual) {
		// Get measurements
		double cost = (double) individual.getOverallCost() / budget;
		double coverage = (double) individual.coverage() / CountryMap.countrySize;
		
		// Logarithmic raise for good coverage or too expensive distributions
		// Graph at http://goo.gl/X0LguK
		// cost = - Math.log(Math.max(0, 1.0 - 0.9 * cost));
		// coverage = - Math.log(Math.max(0, 1.0 - 0.9 * coverage));
		
		// Sigmoid alternative
		// Graph at http://goo.gl/n384EP
		// coverage = 1 / (1 + Math.exp(-8 * coverage + 4));
		
		// Linear relation
		// Graph at http://goo.gl/bS4OZ3
		// double fitness = coverage / (1 + cost);
		
		// Clamp after goal reached
		cost = clamp(cost, 1.0, Integer.MAX_VALUE) - 1.0;
		coverage = clamp(coverage, 0.0, 0.9) / 0.9;
		double fitness = coverage / (1 + cost);
		
		return (float) fitness;
	}
	
	public MastDistribution mutate(MastDistribution individual, int generation) {
		double rate = 1 / Math.pow(1 + (double) generation / 1000, 5);
		rate = clamp(rate, 0.01, 1.0);
		rate = 0.5;
		
		// Get masts as hashmap
		HashMap<Point, TransmitterMast> masts = getMasts(individual);
		
		// Perform random manipulations
		Random random = new Random(System.currentTimeMillis());
		switch (random.nextInt(3)) {
		case 0: changeMastAmount(masts, 0.1 * rate * random()); break;
		case 1: changeMastPositions(masts, 3 - generation * 3 / 1000); break;
		case 2: changeMastTypes(masts, 0.5 * rate * Math.random()); break;
		}
		
		// Create individual from new mast set
		MastDistribution successor = new MastDistribution();
		setMasts(successor, masts);
		return successor;
	}
	
	public MastDistribution crossover(MastDistribution... parentIndividuals) {
		Random random = new Random(System.currentTimeMillis());
		MastDistribution result = new MastDistribution();
		for (int y = 0; y < CountryMap.mapHeight; ++y) {
			for (int x = 0; x < CountryMap.mapWidth; ++x) {
				int choice = random.nextInt(parentIndividuals.length);
				TransmitterMast mast = parentIndividuals[choice].getMastMap()[y][x];
				if (mast != null)
					try {
						result.addMast(mast, x, y);
					} catch (Exception e) {
						e.printStackTrace();
					}
			}
		}
		return result;
	}
	
	public boolean isTerminated() {
		//TODO: implement
		//e.g., stop after first solution (solution!=null) or try to find better ones?
		return false;
	}
	
	private HashMap<Point, TransmitterMast> getMasts(MastDistribution individual) {
		HashMap<Point, TransmitterMast> masts = new HashMap<Point, TransmitterMast>();		
		for (int y = 0; y < individual.getMastMap().length; ++y) {
			TransmitterMast[] row = individual.getMastMap()[y];
			for (int x = 0; x < row.length; ++x) {
				TransmitterMast mast = row[x];
				if (mast != null)
					masts.put(new Point(x, y), mast);
			}
		}
		return masts;
	}

	private void setMasts(MastDistribution individual, HashMap<Point, TransmitterMast> masts) {
		individual.clear();
		for (Point point : masts.keySet()) {
			try {
				individual.addMast(masts.get(point), point.x, point.y);
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
	}
	
	private void addMasts(HashMap<Point, TransmitterMast> masts, int amount) {
	    Random random = new Random(System.currentTimeMillis());
	    
	    // Add given amount of masts
		for (int i = 0; i < amount; ++i) {
			// Decide for mast type
			TransmitterMast mast = randomMast();
			
			// Find valid coordinates
			for (int abort = 0; abort < 10; ++abort) {
				int x = random.nextInt(CountryMap.mapWidth);
				int y = random.nextInt(CountryMap.mapHeight);
				Point point = new Point(x, y);
				if (CountryMap.countryMap[y][x] != 0 && masts.get(point) == null) {
					masts.put(point, mast);
					break;
				}
			};
		}
	}

	private TransmitterMast randomMast() {
		Random random = new Random(System.currentTimeMillis());
		switch (random.nextInt(3)) {
			case 0: return new TransmitterMast(TransmitterMast.MastType.A);
			case 1: return new TransmitterMast(TransmitterMast.MastType.B);
			case 2: return new TransmitterMast(TransmitterMast.MastType.C); 
		}
		return null;
	}
	
	private Point pickRandomMast(HashMap<Point, TransmitterMast> masts) {
		Random random = new Random(System.currentTimeMillis());
		Point[] points = masts.keySet().toArray(new Point[masts.size()]);
		return points[random.nextInt(points.length)];
	}

	private void changeMastAmount(HashMap<Point, TransmitterMast> masts, double rate) {
		// Change amount of masts
		int amount = (int) (masts.size() * (1 + rate));
		amount = Math.max(amount, 1);
		while (amount < masts.size())
			masts.remove(pickRandomMast(masts));
		
		// Fill up with random masts
		addMasts(masts, amount - masts.size());
	}
	
	private void changeMastPositions(HashMap<Point, TransmitterMast> masts, int distance) {
		// For each each point
		@SuppressWarnings("unchecked")
		HashMap<Point, TransmitterMast> from = (HashMap<Point, TransmitterMast>) masts.clone();
		masts.clear();
		for (Point point : from.keySet()) {
			// Try ten times
			for (int i = 0; i < 10; ++i) {
				// Find a near place
				int x = point.x + (int) random(0.5-distance, 0.5+distance);
				int y = point.y + (int) random(0.5-distance, 0.5+distance);
				Point place = new Point(x, y);
				
				// Skip if out of country or already blocked
				if (0 > x || x > CountryMap.mapWidth - 1)
					continue;
				if (0 > y || y > CountryMap.mapHeight - 1)
					continue;
				if (CountryMap.countryMap[y][x] == 0)
					continue;
				if (masts.get(place) != null)
					continue;
	
				// Move mast to new place
				masts.put(place, from.get(point));
				break;
			}
		}
	}
	
	private void changeMastTypes(HashMap<Point, TransmitterMast> masts, double percentage) {
		int changes = (int) (masts.size() * percentage);
		for (int i = 0; i < changes; ++i)
			masts.put(pickRandomMast(masts), randomMast());
	}
	
	private double random() {
		return random(-1.0, 1.0);
	}
	
	private double random(double min, double max) {
		return min + (max - min) * Math.random();
	}
	
	private double clamp(double value) {
		return clamp(value, 0.0, 1.0);	
	}
	
	private double clamp(double value, double min, double max) {
		return Math.min(max, Math.max(min, value));
	}
	
	private int sign(double value) {
		return value < 0.0 ? -1 : value > 0.0 ? 1 : 0;
	}

	public static void main(String[] args) {
		long budget = 100000000;
		float minCoverage = 0.9f;
		// The number of individuals in a generation. You may change this parameter
		// if you need to. However, if you choose to change it - leave a comment why.
		int generationSize = 1000;
		Evolution evolution = new Evolution(generationSize, budget, minCoverage);
		
		//stop after "terminated signal or after n iterations (default: 1000)
		while (evolution.isTerminated() && evolution.getGenerationNumber() < 1000) { 
			evolution.evolve();
    	}
	}
}
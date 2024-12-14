
import { unstable_cache } from 'next/cache';
import { revalidateTag } from 'next/cache';

///////////////////////////////////////////////////////
// type definition section
// Define the type for all possible JSON values
// type JSONPrimitive = string | number | boolean | null; // Primitive JSON types
// type JSONValue = JSONPrimitive | JSONValue[] | { [key: string]: JSONValue }; // Recursive definition for objects and arrays
interface AntibioResMutationItem {
    id: number; // or string, depending on your API response
    mutation: string;
}

// Define the type for the data returned by fetchData
type DataType = AntibioResMutationItem[]; //'any'; // Replace 'any' with your actual data type


// Define an interface for the cache configuration
interface CacheConfig<DataType> {
    fetchData: (url: string) => Promise<DataType>; // Function type for fetching data
    keyParts: string[];                             // Array of strings for key parts
    options: CacheOptions;                          // Options using the CacheOptions interface
}

// Define an interface for the options parameter
interface CacheOptions {
    revalidate?: number; // Optional revalidation interval in seconds
    tags?: string[];     // Optional array of tags for cache invalidation
}

//////////////////////////////////////////////////////

async function fetchData(url: string):  Promise<DataType> {
    // a debug hack, to convince oneself this happens on the server side
    // MYHOSTNAME must be define as the env variable at the place the server is started
    const hostname = process.env.MYHOSTNAME;
    console.log(`**** Greetings of the day. I am running on ${hostname}.`);
    const response = await fetch(url);
    if (!response.ok) {
        throw new Error('Network response was not ok.');
    }
    return response.json();
}


// Create a cached version of fetchData using unstable_cache with the defined interface
const cacheConfig: CacheConfig<DataType> = {
    fetchData,
    keyParts: ['external-data'],
    options: {
        revalidate: 2, // Cache for two secs () or any time you like
        tags: ['external-data'],  // Tag for invalidation
    },
};

const cachedFetchData = unstable_cache(
    cacheConfig.fetchData,
    cacheConfig.keyParts,
    cacheConfig.options
);


export default async function Home() {
    let menuItems: AntibioResMutationItem[] = [];
    try {
        menuItems = await cachedFetchData('http://127.0.0.1:8000/bad_bac_exercise/arm/');
        // menuItems = await cachedFetchData('http://baxter-backend-ctnr:8000/bad_bac_exercise/arm/');
    } catch (err) {
        console.error(err);
    }
    if (menuItems === null) {
        menuItems = []
    }
   // Invalidate cache for external data - here, this is for demo purposes
    revalidateTag('external-data');

    return (
        <div>
            <h1>Menu Items</h1>
             <ul>
                {menuItems.map(item => (
                    <li key={item.id}>{item.mutation}</li>
                ))}
            </ul>
       </div>
    );
}
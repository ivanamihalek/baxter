For production, make sure that `next-config.mjs`
specifies `standalone` output
```javascript
/** @type {import('next').NextConfig} */
const nextConfig = {
    output: "standalone",
};

export default nextConfig;
```

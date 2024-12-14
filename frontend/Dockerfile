# Base image with Node.js
# The FROM instruction initializes a new build stage and sets the base image for subsequent instructions.
FROM node:18-alpine AS build
# Install dependencies only when needed
# Check https://github.com/nodejs/docker-node/tree/b4117f9333da4138b03a546ec926ef50a31506c3#nodealpine
#  to understand why libc6-compat might be needed.
RUN apk add --no-cache libc6-compat

# Install pnpm globally
RUN npm install -g pnpm@9

# The WORKDIR instruction sets the working directory for any
# RUN, CMD, ENTRYPOINT, COPY and ADD instructions that follow.
WORKDIR /app
# Copy and install the dependencies for the project
COPY package.json pnpm-lock.yaml ./
# RUN npm ci
# Copy all other project files to working directory
COPY . .
# Run the next build process and generate the artifacts
RUN pnpm install --frozen-lockfile && pnpm build

# we are using multi stage build process to keep the image size as small as possible
FROM node:18-alpine
# update and install latest dependencies, add dumb-init package
# add a non root user
RUN apk update && apk upgrade && apk add dumb-init && adduser -D nextuser

# set work dir as app
WORKDIR /app
# copy the public folder from the project as this is not included in the build process
# COPY --from=build --chown=nextuser:nextuser /app/public ./public
# copy the standalone folder inside the .next folder generated from the build process
COPY --from=build --chown=nextuser:nextuser /app/.next/standalone ./
# copy the static folder inside the .next folder generated from the build process
COPY --from=build --chown=nextuser:nextuser /app/.next/static ./.next/static
# set non root user
USER nextuser

# expose 3000 on container
EXPOSE 3000

# set app host ,port and node env
ENV HOST=0.0.0.0 PORT=3000 NODE_ENV=production
# start the app with dumb init to spawn the Node.js runtime process
# with signal support
CMD ["dumb-init","node","server.js"]


# Command to run the application
# CMD ["node", "server.js"]
# CMD ["pnpm", "start"]

# using CMD ["node", ".next/standalone/server.js"] focuses on running a lightweight,
# optimized version of your Next.js application, while CMD ["pnpm", "start"] runs
# the application through pnpm, which may include additional overhead but allows for more complex setups.
# The choice between the two depends on your specific needs regarding performance, image size, and application complexity.

# [IM: ] In my hands, the standalone version is both bigger,
# and it cannot find some chunks (so its not working)
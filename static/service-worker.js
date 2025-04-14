const CACHE_NAME = 'seng645-cache-v1';
const urlsToCache = [
  '/',
  '/login',
  '/register',
  '/home',
  '/history',
  '/logout',
  '/stream',
  '/static/styles.css',
  '/static/icons/icon-192x192.png',
  '/static/icons/icon-512x512.png'
];

self.addEventListener('install', event => {
  event.waitUntil(
    caches.open(CACHE_NAME)
      .then(cache => cache.addAll(urlsToCache))
  );
});

self.addEventListener('fetch', event => {
  event.respondWith(
    caches.match(event.request).then(response => response || fetch(event.request))
  );
});
